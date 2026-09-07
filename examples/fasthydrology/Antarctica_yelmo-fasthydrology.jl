## Antarctica counterpart of Greenland_example_yelmo-fasthydro.jl -- same coupling recipe
## (FastHydrology K24/HAB/Shakti <-> Yelmo, both backends), different domain. Read that
## script first; this one is deliberately structured identically (same function names, same
## order) so the two are easy to diff. Domain-specific pieces only: input data (191x191 @
## 32km vs 106x181 @ 16km), the Antarctica_mirror.nml Mirror namelist, and PLOT_DIR (its own
## subdirectory, so this never overwrites the Greenland script's plots).
##
## No forcing (smb/T_srf/Q_geo) is set here, same simplification the Greenland script already
## makes -- this is a coupling-mechanics smoke test, not a physically forced run; see
## tutorial_yelmo_antarctica.jl if you want a properly forced Antarctica setup to build from.

## Preamble ############################################################
# Run from examples/: that is where Project.toml/Manifest.toml live, and where the
# Fortran side resolves its relative input/ paths (input/yelmo_defaults.nml,
# input/yelmo_phys_const.nml) and ice_data/. This script lives one level down,
# so cd to the parent rather than to @__DIR__.
cd(dirname(@__DIR__))
import Pkg; Pkg.activate(".")
#######################################################################

using Revise
using CairoMakie
using Oceananigans.Fields
using Oceananigans.Grids
using Oceananigans.AbstractOperations: Average
using Yelmo
using Yelmo: YelmoMirror                       # Fortran-backed backend: same physics, same nml,
                                                 # same data -- ccall into libyelmo_c_api.so instead
                                                 # of running Julia code directly.
using Yelmo.YelmoPar: YelmoParameters           # YelmoMirror's parameter type (structurally the
                                                 # same nml blocks as YelmoModelParameters, plus
                                                 # p.phys since Mirror keeps constants Fortran-side).
using IceSheetBenchmarks
using Statistics
using FastHydrology
using Shakti

# ── Constants ─────────────────────────────────────────────────────────────────

# CODING_DIR defaults to this machine's layout but is overridable (e.g. on a cluster where
# sibling repos live somewhere else) via the CODING_DIR environment variable.
const CODING_DIR                   = get(ENV, "CODING_DIR", "/Users/taange001/Documents/Coding")
const YELMO_PROJECT                = joinpath(CODING_DIR, "Yelmo.jl")
const ICE_SHEET_BENCHMARKS_PROJECT = joinpath(YELMO_PROJECT, "benchmarks/IceSheetBenchmarks")
const RUN_DIR                      = @__DIR__
# Unlike Greenland's InitMIP-GRL data (committed inside Yelmo.jl/benchmarks/), no equivalent
# official Antarctica benchmark dataset exists there -- this uses the same locally-transferred
# ice_data/Antarctica/ANT-32KM/ files tutorial_yelmo_antarctica.jl uses (see that script's
# header comment: BedMachine topography, RACMO climate, no velocity/age product transferred).
const DATA_DIR                     = joinpath(dirname(RUN_DIR), "ice_data", "Antarctica", "ANT-32KM")
const PLOT_DIR                     = joinpath(dirname(RUN_DIR), "plots", "antarctica")   # own subdirectory,
    # so this never collides with Greenland_example_yelmo-fasthydro.jl's plots/*.png

# Reuses run01/Greenland.nml's physics/solver settings verbatim -- same simplification
# tutorial_yelmo_antarctica.jl makes (see its header comment): every knob it sets is
# domain-agnostic Yelmo config, so this only has to teach the *data* side of a new domain.
const YELMO_NML          = joinpath(dirname(RUN_DIR), "run01", "Greenland.nml")   # lives inside this repo, so always relative to it
const YELMO_REGIONS_FILE = joinpath(DATA_DIR, "ANT-32KM_REGIONS.nc")
const YELMO_BASINS_FILE  = joinpath(DATA_DIR, "ANT-32KM_BASINS-nasa.nc")
const YELMO_TOPO_FILE    = joinpath(DATA_DIR, "ANT-32KM_TOPO-BedMachine.nc")

# Yelmo works in year-based units (velocities in m/yr, ATT in Pa^-3 yr^-1) while FastHydrology
# (and Shakti) work in SI/second-based units, so every field crossing the coupling needs
# perYear2perSecond. rho_i must match KazmierczakHydroModel's/Shakti's own rho_i default, since
# it's used to convert Yelmo's ice-equivalent bmb_grnd into a mass melt rate for K24.
const RHO_I = 917.0

# Outer (Yelmo) timestep and total run length, shared by every coupling branch. K24/HAB are
# steady-state (no internal dt of their own -- they just re-solve from whatever Yelmo state is
# pushed in each step), so DT_YR only really matters for Shakti, which is genuinely time-evolving:
# its own internal dt (see build_hydrology_sim_Shakti) is set equal to DT_YR (in seconds) rather
# than sub-cycled at Shakti's own finer native timestep, so one Shakti step covers exactly one
# Yelmo step. This is deliberately not physically resolved for Shakti's fast subglacial-flow
# timescale -- kept cheap (few steps) purely to check the push/step/pull cycle runs end to end and
# every field is properly refreshed, not to produce a physically converged Shakti run.
const DT_YR       = 2.0
const TIME_END_YR = 4.0

# ── Backend selection ────────────────────────────────────────────────────────
# "yelmo" (default) — pure-Julia YelmoModel, as this script has always run.
# "mirror"           — YelmoMirror, the Fortran-backed backend (ccall into libyelmo_c_api.so).
#                       Set FASTHYDRO_BACKEND=mirror to use this instead.
#
# Every coupling function below (Yelmo_to_FastHydrology_*!, FastHydrology_to_Yelmo_*!,
# couple_step!, step!, run!) is reused unchanged across both backends -- they only touch
# interior(...) fields, which both backends expose identically. What differs is model
# construction (build_yelmo vs build_yelmo_mirror) and how the "N_eff is set externally" flag
# is expressed: yneff.method = -1 in Julia's YelmoModelParameters vs. hyd.is_external = true in
# Fortran's &yhyd (see yelmo/src/yelmo_dynamics.f90::calc_ydyn_neff's hyd%par%is_external
# bypass, added specifically so an external host can push N_eff into YelmoMirror the same way).
const BACKEND = lowercase(get(ENV, "FASTHYDRO_BACKEND", "yelmo"))
BACKEND in ("yelmo", "mirror") ||
    error("FASTHYDRO_BACKEND must be 'yelmo' or 'mirror', got '$BACKEND'")

# Mirror-only: own nml (not YELMO_NML domain-overridden) since Antarctica's real datasets use
# different filename templates (BedMachine/RACMO, not M17/MARv3.11) that the nml's
# ice_data/{domain}/{grid_name}/... paths encode directly, and its &yelmo_data pd_*_load
# flags are turned off (no velocity/age product transferred for this domain) -- see that
# file's own header comment. Shared with tutorial_yelmo_antarctica.jl, not duplicated here.
const YELMO_NML_MIRROR = joinpath(dirname(RUN_DIR), "tutorial", "Antarctica_mirror.nml")

# Mirror-only workaround: YelmoMirror's Fortran-side yelmo_init cannot safely be called a
# second time in the same process when building from a file-based grid (grid=nothing, as
# build_yelmo_mirror does below) -- observed to crash in ydyn_alloc ("allocatable array is
# already allocated") on the second of four sequential branch builds. This is a pre-existing
# Fortran-side issue unrelated to the hydrology coupling this script tests (crashes allocating
# plain velocity-diagnostic arrays, nothing hydrology- or N_eff-specific), not root-caused here
# -- out of scope for wiring up the coupling. Workaround: run one branch per process. When set,
# FASTHYDRO_BRANCH restricts main()'s loop to that single 1-based branch index; unset (the
# "yelmo" backend's normal mode) runs all four branches in one process, as before.
const BRANCH_SEL = get(ENV, "FASTHYDRO_BRANCH", "")

# ── Helpers ───────────────────────────────────────────────────────────────────

"""Replace field `f` in struct `s` with value `v`, returning a new instance."""
function _override_field(s, f::Symbol, v)
    vals = (n === f ? v : getfield(s, n) for n in fieldnames(typeof(s)))
    return typeof(s)(vals...)
end

"""Compute a lazy `BinaryOperation` and return the interior 2-D slice."""
function _compute_interior(op)
    result = Field(op)
    compute!(result)
    return interior(result, :, :, 1)
end

# ── Yelmo setup ───────────────────────────────────────────────────────────────

"""`external_neff = true` overrides `yneff.method` to -1 (see `setup_yelmo` for the full
rationale) -- only appropriate for a branch where hydrology coupling will actually push N_eff
every step. The "No coupling" branch must pass `external_neff = false` and keep the shared nml's
own method = 3 (van Pelt & Bueler till closure): nothing would ever push N_eff for it, so under
method = -1 it would stay at its Field-allocation default (zero) for the whole run -- a
frictionless bed."""
function build_yelmo_parameters(; external_neff::Bool)
    p = YelmoModelParameters(YELMO_NML, "Antarctica")

    # YELMO_NML (shared with the Greenland script) has domain = "Greenland" baked into its
    # &yelmo block; override it explicitly since nothing else here re-derives it from the
    # "Antarctica" name label above (matches tutorial_yelmo_antarctica.jl's build_yelmo_parameters).
    yelmo = _override_field(p.yelmo, :domain, "Antarctica")
    p     = _override_field(p, :yelmo, yelmo)

    init_topo = _override_field(p.yelmo_init_topo, :init_topo_path, YELMO_TOPO_FILE)

    masks = _override_field(p.yelmo_masks, :regions_path, YELMO_REGIONS_FILE)
    masks = _override_field(masks,         :basins_path,  YELMO_BASINS_FILE)

    p = _override_field(p, :yelmo_init_topo, init_topo)
    p = _override_field(p, :yelmo_masks,     masks)

    if external_neff
        # run01/Greenland.nml sets yneff.method = 3 (van Pelt & Bueler till closure), which makes
        # calc_ydyn_neff! *recompute* dyn.N_eff from thrm.H_w every step -- silently discarding
        # whatever N_eff the coupling below just pushed in from K24/HAB/Shakti before Yelmo.step!
        # ever uses it (H_w still feeds through, but as an input to Yelmo's own till closure, not
        # as the hydrology model's own N). method = -1 is the one setting that makes
        # calc_ydyn_neff! a no-op ("N_eff set externally -- leave alone", see
        # Yelmo.jl/src/dyn/neff.jl), so the pushed N_eff is what the velocity solve actually sees.
        yneff = _override_field(p.yneff, :method, -1)
        p     = _override_field(p, :yneff, yneff)
    end

    return p
end

function build_yelmo(; external_neff::Bool)
    mkpath(RUN_DIR)

    for (path, label) in [
        (YELMO_NML,          "Yelmo namelist"),
        (YELMO_REGIONS_FILE, "regions file"),
        (YELMO_BASINS_FILE,  "basins file"),
        (YELMO_TOPO_FILE,    "topography file"),
    ]
        isfile(path) || error("Yelmo $label not found: $path")
    end

    p         = build_yelmo_parameters(; external_neff)
    benchmark = InitMIPGRLBenchmark(YELMO_REGIONS_FILE)
    y         = YelmoModel(benchmark, 0.0; p, boundaries=:bounded, rundir=RUN_DIR)

    init_topo_load!(y; grad_lim_zb=p.ytopo.grad_lim_zb)
    init_masks!(y)
    fill!(interior(y.bnd.H_sed), 100.0)

    return p, y
end

"""Mirror counterpart of `build_yelmo_parameters`. `external_neff = true` overrides
`hyd.is_external` to `true` -- the &yhyd equivalent of the Julia backend's `yneff.method = -1`
(see `build_yelmo_parameters`'s docstring for the full rationale, and
`yelmo/src/yelmo_dynamics.f90::calc_ydyn_neff` for the Fortran-side bypass this flag controls).
Same contract as the Julia backend: only appropriate for a branch that will actually push
N_eff every step. "No coupling" must pass `external_neff = false` and keep the shared nml's own
&yhyd default (`is_external = false`), so `calc_ydyn_neff` keeps recomputing `dyn%now%N_eff`
from `hyd%now%N` (FastHydrology's till closure) every step -- nothing ever pushes N_eff for
that branch, same reasoning as `yneff.method = 3` on the Julia side."""
function build_yelmo_parameters_mirror(; external_neff::Bool)
    p = YelmoParameters(YELMO_NML_MIRROR, "Antarctica")

    if external_neff
        hyd = _override_field(p.hyd, :is_external, true)
        p   = _override_field(p, :hyd, hyd)
    end

    return p
end

"""Mirror counterpart of `build_yelmo`. Unlike the Julia backend, Mirror's Fortran-side
`yelmo_init` reads topography/masks/regions itself from the nml's own
`ice_data/{domain}/{grid_name}/...` path templates, resolved relative to the process's current
working directory (`RUN_DIR`, since the preamble `cd`s there) -- see
`examples/ice_data/Antarctica/ANT-32KM/` for the symlinks this relies on (same data as
DATA_DIR, under the exact filenames the nml expects). No `init_topo_load!`/`init_masks!`
equivalent is needed here; Fortran does that itself during `yelmo_init`."""
function build_yelmo_mirror(; external_neff::Bool)
    mkpath(RUN_DIR)
    isfile(YELMO_NML_MIRROR) || error("Yelmo Mirror namelist not found: $YELMO_NML_MIRROR")

    p = build_yelmo_parameters_mirror(; external_neff)
    y = YelmoMirror(p, 0.0; alias="fasthydro_mirror_ant", rundir=RUN_DIR, overwrite=true)
    fill!(interior(y.bnd.H_sed), 100.0)

    return p, y
end

# ── Coupling type hierarchy ───────────────────────────────────────────────────

abstract type HydrologyCoupling end

"""Fully coupled: Yelmo ↔ FastHydrology exchange each timestep. `S` is whichever FastHydrology
simulation type (`SteadyStateSimulation{KazmierczakHydroModel}`, `SteadyStateSimulation{HABHydroModel}`,
or `TimeSimulation{ShaktiHydroModel}`) `build_hydrology_sim_*` below produced -- `couple_step!`
dispatches on `sim.model`'s type to run the right coupling for whichever one is active."""
struct CoupledHydrology{S} <: HydrologyCoupling
    sim::S
end

"""No hydrology: Yelmo advances alone."""
struct NoCoupling <: HydrologyCoupling end

# ── K24 (Kazmierczak et al. 2024): build + coupling functions ─────────────────

"""Build a FastHydrology SteadyStateSimulation wrapping KazmierczakHydroModel, from `yelmo`'s current fields."""
function build_hydrology_sim_K24(yelmo)

    T = Float32

    Nx, Ny = yelmo.g.Nx, yelmo.g.Ny
    xlims  = (0, T(yelmo.g.Δxᶜᵃᵃ * Nx))
    ylims  = (0, T(yelmo.g.Δyᵃᶜᵃ * Ny))

    mask    = _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd) .> 0
    h       = interior(yelmo.tpo.H_ice,  :, :, 1)   # ice thickness
    b       = interior(yelmo.bnd.z_bed,  :, :, 1)   # bedrock elevation
    abs_v_b = perYear2perSecond.(interior(yelmo.dyn.uxy_b, :, :, 1))              # basal speed
    A_visc  = perYear2perSecond.(mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1])  # depth-averaged rate factor
    # bmb_grnd is ice-equivalent [m/yr], positive = accretion, negative = melt; mdot wants a
    # positive mass melt rate [kg/m^2/s], hence the sign flip and the rho_i scaling.
    mdot    = perYear2perSecond.(-interior(yelmo.thrm.bmb_grnd, :, :, 1) .* RHO_I)
    kappa   = zeros(T, Nx, Ny)   # bed hardness (0: hard, 1: soft)

    longcoupwater = 0.0   # smoothing of geometric-potential gradients
    fill_iters    = 10    # iterations to fill local minima in potential field

    grid  = OGRectHydroGrid(Nx, Ny, xlims, ylims; T = T)
    model = KazmierczakHydroModel(grid, kappa, abs_v_b, A_visc, mdot; longcoupwater=longcoupwater, fill_iters=fill_iters)
    state = HydroState(grid, mask, h, b)
    return SteadyStateSimulation(model, grid, state)
end

"""Copy Yelmo fields → FastHydrology state/model (K24)."""
function Yelmo_to_FastHydrology_K24!(sim, yelmo)
    sim.state.mask    .= _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd) .> 0
    sim.state.h       .= interior(yelmo.tpo.H_ice, :, :, 1)
    sim.state.b       .= interior(yelmo.bnd.z_bed, :, :, 1)
    sim.model.abs_v_b .= perYear2perSecond.(interior(yelmo.dyn.uxy_b, :, :, 1))
    sim.model.A_visc  .= perYear2perSecond.(mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1])
    sim.model.mdot    .= perYear2perSecond.(-interior(yelmo.thrm.bmb_grnd, :, :, 1) .* RHO_I)
    fill!(sim.model.kappa, 0)
end

"""Copy FastHydrology outputs → Yelmo boundary fields (K24)."""
function FastHydrology_to_Yelmo_K24!(yelmo, sim)
    yelmo.dyn.N_eff .= sim.state.N   # effective pressure
    yelmo.thrm.H_w  .= sim.state.W   # subglacial water thickness
end

# ── HAB (height above buoyancy): build + coupling functions ───────────────────

"""Build a FastHydrology SteadyStateSimulation wrapping HABHydroModel, from `yelmo`'s current fields.
HAB needs no model-specific input fields (no abs_v_b/A_visc/mdot/kappa), only the shared
mask/h/b state."""
function build_hydrology_sim_HAB(yelmo)

    T = Float32
    Nx, Ny = yelmo.g.Nx, yelmo.g.Ny
    xlims  = (0, T(yelmo.g.Δxᶜᵃᵃ * Nx))
    ylims  = (0, T(yelmo.g.Δyᵃᶜᵃ * Ny))

    mask = _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd) .> 0
    h    = interior(yelmo.tpo.H_ice, :, :, 1)
    b    = interior(yelmo.bnd.z_bed, :, :, 1)

    grid  = OGRectHydroGrid(Nx, Ny, xlims, ylims; T = T)
    model = HABHydroModel(grid)
    state = HydroState(grid, mask, h, b)
    return SteadyStateSimulation(model, grid, state)
end

"""Copy Yelmo fields → FastHydrology state (HAB)."""
function Yelmo_to_FastHydrology_HAB!(sim, yelmo)
    sim.state.mask .= _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd) .> 0
    sim.state.h    .= interior(yelmo.tpo.H_ice, :, :, 1)
    sim.state.b    .= interior(yelmo.bnd.z_bed, :, :, 1)
end

"""Copy FastHydrology outputs → Yelmo boundary fields (HAB has no time-evolving water layer, only N)."""
function FastHydrology_to_Yelmo_HAB!(yelmo, sim)
    yelmo.dyn.N_eff .= sim.state.N
end

# ── Shakti: build + coupling functions ────────────────────────────────────────

"""Build a FastHydrology TimeSimulation wrapping a Shakti simulation, seeded from `yelmo`'s current
fields. Modeling choices with no direct Yelmo equivalent (flagged inline): Shakti's mask
categories, the initial gap height, and no moulin/point-source coupling (melt is generated
internally from geothermal flux + frictional heating, matching Shakti's own real-glacier example).

`dt_yr` (Yelmo's own outer timestep, in years) becomes Shakti's internal dt directly (converted to
seconds) rather than a finer native Shakti timestep that couple_step! would sub-cycle -- one Shakti
step per Yelmo step. Not physically resolved for Shakti's own (much faster) timescale, but this
coupling is a smoke test of the push/step/pull cycle, not a physically converged Shakti run; see
the DT_YR/TIME_END_YR comment above."""
function build_hydrology_sim_Shakti(yelmo, dt_yr)

    T = Float64   # Shakti's own solvers are written for Float64

    Nx, Ny = yelmo.g.Nx, yelmo.g.Ny
    dx, dy = yelmo.g.Δxᶜᵃᵃ, yelmo.g.Δyᵃᶜᵃ

    grid = Shakti.Grid(Nx, Ny, T((Nx - 1) * dx), T((Ny - 1) * dy))

    # Shakti's GROUNDED is the same physical concept as FastHydrology's HydroState.mask == 1
    # (grounded ice present) -- just one of four categories here (vs a boolean there), since
    # Shakti also needs to distinguish ocean/land/other-basin among the non-grounded cells that
    # K24/HAB never look at. Built independently (not derived from a HydroState), and static for
    # the whole run -- set only here, not per step.
    grounded_frac = _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd)
    zb0 = interior(yelmo.bnd.z_bed, :, :, 1)
    mask = fill(Shakti.OTHER_BASIN, Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        if grounded_frac[i, j] >= 0.5
            mask[i, j] = Shakti.GROUNDED
        elseif zb0[i, j] < 0
            mask[i, j] = Shakti.OCEAN
        else
            mask[i, j] = Shakti.LAND
        end
    end
    mask[[1, Nx], :] .= Shakti.OTHER_BASIN   # close off the domain edge (closed-boundary convention)
    mask[:, [1, Ny]] .= Shakti.OTHER_BASIN

    zb     = zb0
    H      = interior(yelmo.tpo.H_ice, :, :, 1)
    # Shakti derives H = zs - zb - b internally (compute_H!), so zs must be reconstructed as
    # zb + H_ice here rather than read from yelmo.tpo.z_srf directly: z_srf is referenced to sea
    # level (0 over open ocean), not to zb + H_ice, so feeding it in directly would make Shakti
    # read open-ocean bathymetry as if it were ice thickness (matches Shakti's own real-glacier
    # example, Helheim.jl, which builds zs the same way).
    zs     = zb .+ H
    A_visc = perYear2perSecond.(mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1])
    ub_x   = perYear2perSecond.(interior(yelmo.dyn.ux_b, :, :, 1))
    ub_y   = perYear2perSecond.(interior(yelmo.dyn.uy_b, :, :, 1))
    taub_x = interior(yelmo.dyn.taub_acx, :, :, 1)   # already in Pa, no time unit to convert
    taub_y = interior(yelmo.dyn.taub_acy, :, :, 1)
    G      = interior(yelmo.bnd.Q_geo, :, :, 1) .* 1e-3   # mW/m^2 -> W/m^2 - WARNING: this field might contain a -9999.0 fill sentinel
    gap0   = fill(1e-3, Nx, Ny)   # initial gap height guess; no Yelmo equivalent
    ieb    = zeros(Nx, Ny)        # no explicit moulin coupling (see docstring above)

    p  = Shakti.ModelParameters(rho_i = RHO_I)
    mi = Shakti.ConstantMeltInput()
    sl = Shakti.PrescribedSlidingLaw()   # taub is supplied directly from Yelmo's own dynamics solve, not solved by Shakti

    state = Shakti.State(grid)
    Shakti.set_initial_conditions!(state, grid, p, sl, mask, A_visc, zb, zs, gap0, G, ub_x, ub_y, ieb, taub_x, taub_y)

    ls = Shakti.CholeskyDirectSolver(grid)
    ps = Shakti.PicardSolver(500, 1e-6, ls, grid)

    # tsteps is inert here: it only matters for Shakti.run!'s own loop/observer bookkeeping, and
    # this coupling drives Shakti with FastHydrology.step! (see couple_step! below), not run!.
    dt_seconds = dt_yr * FastHydrology.SECONDS_PER_YEAR
    shakti_sim = Shakti.Simulation(grid, state, 1, dt_seconds, p, "implicit", String[], mi, sl; ps = ps)

    model = ShaktiHydroModel(shakti_sim)
    return TimeSimulation(model)
end

"""Push Yelmo's current ice geometry, rheology, and basal stress/velocity into the wrapped Shakti
simulation. Shakti's mask is static per run (set only in build_hydrology_sim_Shakti), so it is not
touched here. Unlike `set_initial_conditions!`, Shakti's own per-timestep `step!` does *not*
re-derive H/po/abs_ub from zs/zb/ub_x/ub_y on its own (only `set_initial_conditions!` and
`compute_beta!` -- called internally each step -- touch those), so this has to call the matching
`compute_*!` refreshers explicitly after pushing the raw fields, or the coupling would silently
keep using the ice geometry/velocity from `build_hydrology_sim_Shakti`'s initial call forever."""
function Yelmo_to_FastHydrology_Shakti!(shakti_sim, yelmo)
    s = shakti_sim.state

    s.zb .= interior(yelmo.bnd.z_bed, :, :, 1)
    H_ice = interior(yelmo.tpo.H_ice, :, :, 1)
    # Reconstructed as zb + H_ice, not read from yelmo.tpo.z_srf -- see build_hydrology_sim_Shakti's
    # comment: z_srf is sea-level-referenced over open ocean, not zb + H_ice, and Shakti's own
    # H = zs - zb - b derivation would otherwise read ocean depth as ice thickness there.
    s.zs .= s.zb .+ H_ice

    s.A_visc .= perYear2perSecond.(mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1])
    s.ub_x   .= perYear2perSecond.(interior(yelmo.dyn.ux_b, :, :, 1))
    s.ub_y   .= perYear2perSecond.(interior(yelmo.dyn.uy_b, :, :, 1))
    s.taub_x .= interior(yelmo.dyn.taub_acx, :, :, 1)
    s.taub_y .= interior(yelmo.dyn.taub_acy, :, :, 1)
    s.G      .= interior(yelmo.bnd.Q_geo, :, :, 1) .* 1e-3

    Shakti.apply_mask_to_sliding!(s)   # re-zero ub_x/ub_y on any face touching an OTHER_BASIN cell
    Shakti.compute_abs_ub!(s)          # |u_b| from the new ub_x/ub_y
    Shakti.compute_H!(s)               # ice thickness from the new zs/zb/b
    Shakti.compute_po!(s, shakti_sim.p) # overburden pressure from the new H
end

"""Copy Shakti's effective pressure / gap height back into Yelmo boundary fields. Shakti's state
arrays are plain (Nx, Ny) arrays (unlike K24/HAB's Oceananigans fields), hence writing through
`interior(...)` on the Yelmo side."""
function FastHydrology_to_Yelmo_Shakti!(yelmo, shakti_sim)
    interior(yelmo.dyn.N_eff, :, :, 1) .= shakti_sim.state.N
    interior(yelmo.thrm.H_w,  :, :, 1) .= shakti_sim.state.b   # gap height stands in for water-layer thickness
end

# ── Coupling dispatch (per active FastHydrology model) ────────────────────────

"""Advance the hydrology model by one Yelmo timestep and push results back, dispatched on
whichever FastHydrology model is active in `sim`."""
function couple_step!(::KazmierczakHydroModel, sim, yelmo, dt)
    Yelmo_to_FastHydrology_K24!(sim, yelmo)
    FastHydrology.run!(sim)
    FastHydrology_to_Yelmo_K24!(yelmo, sim)
end

function couple_step!(::HABHydroModel, sim, yelmo, dt)
    Yelmo_to_FastHydrology_HAB!(sim, yelmo)
    FastHydrology.run!(sim)
    FastHydrology_to_Yelmo_HAB!(yelmo, sim)
end

function couple_step!(model::ShaktiHydroModel, sim, yelmo, dt)
    shakti_sim = model.sim
    # shakti_sim.dt was set to dt_yr * SECONDS_PER_YEAR at construction (build_hydrology_sim_Shakti)
    # -- if the caller ever runs this coupling with a different outer dt than it was built with,
    # one Shakti step would silently cover the wrong span of Yelmo time. Guard rather than let that
    # drift silently.
    dt_seconds = dt * FastHydrology.SECONDS_PER_YEAR
    isapprox(dt_seconds, shakti_sim.dt; rtol=1e-8) || error(
        "couple_step! (Shakti): outer dt = $dt yr ($dt_seconds s) does not match " *
        "shakti_sim.dt = $(shakti_sim.dt) s -- rebuild the Shakti sim with the dt actually used by run!.")

    Yelmo_to_FastHydrology_Shakti!(shakti_sim, yelmo)
    FastHydrology.step!(sim)   # one Shakti step == one Yelmo step (dt matched, see above)
    FastHydrology_to_Yelmo_Shakti!(yelmo, shakti_sim)
end

# ── Time stepping (dispatched on coupling type) ───────────────────────────────

"""Single timestep with active hydrology coupling."""
function step!(coupling::CoupledHydrology, yelmo, t, dt)
    couple_step!(coupling.sim.model, coupling.sim, yelmo, dt)
    Yelmo.step!(yelmo, dt)
end

"""Single timestep, hydrology disabled — Yelmo advances alone."""
function step!(::NoCoupling, yelmo, t, dt)
    Yelmo.step!(yelmo, dt)
end

"""Run the simulation to `time_end` using whatever coupling strategy is provided."""
function run!(coupling::HydrologyCoupling, yelmo; dt=1.0, time_end=3.0)
    for t in dt:dt:time_end
        step!(coupling, yelmo, t, dt)
        # write_output!(yelmo_out, yelmo)
        println("t = $t yr")
    end
end

"""Build a fresh, initialized Yelmo model (same recipe every time, so runs start from the same t=0
state). `external_neff` must match whatever the coupling built from the returned model will do:
`true` for a branch that pushes hydrology-derived N_eff every step (K24/HAB/Shakti), `false` for
"No coupling" (see `build_yelmo_parameters`/`build_yelmo_parameters_mirror`).

Dispatches on the global `BACKEND` ("yelmo" or "mirror"). The N_eff bootstrap seed below applies
identically to both: `init_state!` (called at the end, for either backend) runs an internal
predictor/corrector spin-up -- including one velocity solve -- before the main time loop has
ever called `couple_step!` to push a real hydrology-derived N_eff. For "mirror", `init_state!`
also does `yelmo_sync!` first thing, pushing whatever's currently in `yelmo.dyn.N_eff` (Julia
buffer) into Fortran before that spin-up runs -- so seeding it here works the same way for both
backends, just via a different underlying push mechanism."""
function setup_yelmo(; external_neff::Bool)
    if BACKEND == "mirror"
        p, yelmo = build_yelmo_mirror(; external_neff)
        rho_ice, g = p.phys.rho_ice, p.phys.g
    else
        p, yelmo = build_yelmo(; external_neff)
        yelmo.bnd.H_sed .= 100.0
        rho_ice, g = yelmo.c.rho_ice, yelmo.c.g
    end

    if external_neff
        # With N_eff set externally (yneff.method = -1 / hyd.is_external = true for
        # external_neff branches), calc_ydyn_neff! is a no-op, so without this N_eff would still
        # be sitting at its Field-allocation default (zero) for that first solve -- a
        # frictionless bed, which is what made the SSA solver diverge before this was added.
        # Seed it with an overburden estimate (same formula as yneff.method == 1) as a
        # physically reasonable placeholder; the first real couple_step! call overwrites it
        # before the first actual *time* step.
        H_ice = interior(yelmo.tpo.H_ice, :, :, 1)
        interior(yelmo.dyn.N_eff, :, :, 1) .= rho_ice * g .* H_ice
    end

    init_state!(yelmo, 0.0; thrm_method="robin-cold")
    return yelmo
end

# ── Plotting (dispatched on simulation type) ──────────────────────────────────

"""Visualize effective pressure N, dispatched on whichever FastHydrology simulation type is active."""
plot_N(sim::SteadyStateSimulation, title; kwargs...) =
    visualize_field(sim.state.N; plot_title = "Effective pressure: " * title, kwargs...)

function plot_N(sim::TimeSimulation{<:ShaktiHydroModel}, title; kwargs...)
    shakti_sim = sim.model.sim
    visualize_field(shakti_sim.grid.x, shakti_sim.grid.y, shakti_sim.state.N; plot_title = "Effective pressure: " * title, kwargs...)
end

function main()

    # ── Run ───────────────────────────────────────────────────────────────────────

    mkpath(PLOT_DIR)

    @info "Running FastHydrology coupling smoke test..." backend=BACKEND dt_yr=DT_YR time_end_yr=TIME_END_YR

    plot_title_list    = ["No coupling", "K24 coupling", "HAB coupling", "Shakti coupling"]
    # external_neff must match each branch's coupling: false for "No coupling" (nothing pushes
    # N_eff, so Yelmo's own till closure must stay active), true for the three hydrology-coupled
    # branches (see build_yelmo_parameters / setup_yelmo).
    external_neff_list = [false, true, true, true]
    make_coupling_list = [
        yelmo -> NoCoupling(),
        yelmo -> CoupledHydrology(build_hydrology_sim_K24(yelmo)),
        yelmo -> CoupledHydrology(build_hydrology_sim_HAB(yelmo)),
        yelmo -> CoupledHydrology(build_hydrology_sim_Shakti(yelmo, DT_YR)),
    ]

    # See BRANCH_SEL's definition above: unset runs all four branches in this one process
    # (the "yelmo" backend's normal mode); set (typically by a per-branch sbatch loop for the
    # "mirror" backend) restricts this run to just that one branch.
    branch_indices = isempty(BRANCH_SEL) ? eachindex(make_coupling_list) : [parse(Int, BRANCH_SEL)]

    for idx in branch_indices
        make_coupling = make_coupling_list[idx]

        # Fresh Yelmo model per branch, so all four start from the same t=0 state --
        # otherwise each branch would continue from wherever the previous one left off.
        yelmo    = setup_yelmo(; external_neff = external_neff_list[idx])
        coupling = make_coupling(yelmo)

        # Write to file
        # yelmo_output_groups = [:tpo, :dyn, :thrm, :mat, :bnd]
        # yelmo_out = init_output(yelmo, joinpath(yelmo.rundir, "yelmo.nc"); selection=OutputSelection(groups=yelmo_output_groups))
        # write_output!(yelmo_out, yelmo)

        # Run the total time simulation
        @time run!(coupling, yelmo; dt=DT_YR, time_end=TIME_END_YR)

        # close(yelmo_out)

        # Visualize fields -- headless-safe: saved to PLOT_DIR rather than displayed, since a
        # cluster run has no display to pop a window on. "_mirror" suffix on the mirror backend
        # keeps its plots alongside (not overwriting) the default "yelmo" backend's, so the two
        # can be compared directly.
        slug = lowercase(replace(plot_title_list[idx], " " => "_")) * (BACKEND == "mirror" ? "_mirror" : "")
        if coupling isa CoupledHydrology
            plot_N(coupling.sim, plot_title_list[idx]; display_flag = false, savefig_path = joinpath(PLOT_DIR, "N_$(slug).png"))
        end
        visualize_field(yelmo.dyn.uxy_s;
            plot_title    = "Surface horizontal velocity magnitude: " * plot_title_list[idx],
            display_flag  = false,
            savefig_path  = joinpath(PLOT_DIR, "uxy_s_$(slug).png"),
        )

    end

end

main()
