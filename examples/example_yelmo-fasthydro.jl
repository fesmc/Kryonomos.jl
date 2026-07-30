## Preamble ############################################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#######################################################################

using Revise
using CairoMakie
using Oceananigans.Fields
using Oceananigans.Grids
using Oceananigans.AbstractOperations: Average
using Yelmo
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
const DATA_DIR                     = joinpath(YELMO_PROJECT, "benchmarks/initmip-grl/data/GRL-16KM")
const PLOT_DIR                     = joinpath(RUN_DIR, "plots")

const YELMO_NML          = joinpath(RUN_DIR, "run01", "Greenland.nml")   # lives inside this repo, so always relative to it
const YELMO_REGIONS_FILE = joinpath(DATA_DIR, "GRL-16KM_REGIONS.nc")
const YELMO_BASINS_FILE  = joinpath(DATA_DIR, "GRL-16KM_BASINS-nasa.nc")
const YELMO_TOPO_FILE    = joinpath(DATA_DIR, "GRL-16KM_TOPO-M17-v5.nc")

# Yelmo works in year-based units (velocities in m/yr, ATT in Pa^-3 yr^-1) while FastHydrology
# (and Shakti) work in SI/second-based units, so every field crossing the coupling needs
# perYear2perSecond. rho_i must match KazmierczakHydroModel's/Shakti's own rho_i default, since
# it's used to convert Yelmo's ice-equivalent bmb_grnd into a mass melt rate for K24.
const RHO_I = 917.0

# Shakti's own internal timestep [s] -- much finer than Yelmo's (years), so couple_step! below
# sub-cycles it to cover one Yelmo timestep. 1800 s (30 min) matches Shakti's own real-glacier
# (Helheim) example.
const SHAKTI_DT = 1800.0

# Cap on Shakti sub-steps per Yelmo timestep. NOT physically meaningful -- a real coupled run
# needs enough sub-steps to actually span the elapsed Yelmo dt (years -> ~tens of thousands of
# 30-min Shakti steps), which is too expensive for a smoke test. This cap exists purely so a
# quick "does the whole push/step/pull cycle run end to end" check finishes fast; remove it (or
# raise it a lot) for a physically meaningful coupled run.
const SHAKTI_MAX_SUBSTEPS_SMOKE_TEST = 10

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

function build_yelmo_parameters()
    p = YelmoModelParameters(YELMO_NML, "Greenland")

    init_topo = _override_field(p.yelmo_init_topo, :init_topo_path, YELMO_TOPO_FILE)

    masks = _override_field(p.yelmo_masks, :regions_path, YELMO_REGIONS_FILE)
    masks = _override_field(masks,         :basins_path,  YELMO_BASINS_FILE)

    p = _override_field(p, :yelmo_init_topo, init_topo)
    p = _override_field(p, :yelmo_masks,     masks)

    return p
end

function build_yelmo()
    mkpath(RUN_DIR)

    for (path, label) in [
        (YELMO_NML,          "Yelmo namelist"),
        (YELMO_REGIONS_FILE, "regions file"),
        (YELMO_BASINS_FILE,  "basins file"),
        (YELMO_TOPO_FILE,    "topography file"),
    ]
        isfile(path) || error("Yelmo $label not found: $path")
    end

    p         = build_yelmo_parameters()
    benchmark = InitMIPGRLBenchmark(YELMO_REGIONS_FILE)
    y         = YelmoModel(benchmark, 0.0; p, boundaries=:bounded, rundir=RUN_DIR)

    init_topo_load!(y; grad_lim_zb=p.ytopo.grad_lim_zb)
    init_masks!(y)
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
internally from geothermal flux + frictional heating, matching Shakti's own real-glacier example)."""
function build_hydrology_sim_Shakti(yelmo)

    T = Float64   # Shakti's own solvers are written for Float64

    Nx, Ny = yelmo.g.Nx, yelmo.g.Ny
    dx, dy = yelmo.g.Δxᶜᵃᵃ, yelmo.g.Δyᵃᶜᵃ

    grid = Shakti.Grid(Nx, Ny, T((Nx - 1) * dx), T((Ny - 1) * dy))

    # Shakti's own mask (0/1/2/3 = grounded/ocean/land/other-basin) is a separate concept from
    # FastHydrology's HydroState.mask, and static for the whole run -- set only here, not per step.
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
    zs     = interior(yelmo.tpo.z_srf, :, :, 1)
    H      = interior(yelmo.tpo.H_ice, :, :, 1)
    A_visc = perYear2perSecond.(mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1])
    ub_x   = perYear2perSecond.(interior(yelmo.dyn.ux_b, :, :, 1))
    ub_y   = perYear2perSecond.(interior(yelmo.dyn.uy_b, :, :, 1))
    taub_x = interior(yelmo.dyn.taub_acx, :, :, 1)   # already in Pa, no time unit to convert
    taub_y = interior(yelmo.dyn.taub_acy, :, :, 1)
    G      = interior(yelmo.bnd.Q_geo, :, :, 1) .* 1e-3   # mW/m^2 -> W/m^2
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
    shakti_sim = Shakti.Simulation(grid, state, 1, SHAKTI_DT, p, "implicit", String[], mi, sl; ps = ps)

    model = ShaktiHydroModel(shakti_sim)
    return TimeSimulation(model)
end

"""Push Yelmo's current ice geometry, rheology, and basal stress/velocity into the wrapped Shakti
simulation. Shakti's mask is static per run (set only in build_hydrology_sim_Shakti), so it is not
touched here."""
function Yelmo_to_FastHydrology_Shakti!(shakti_sim, yelmo)
    shakti_sim.state.zb     .= interior(yelmo.bnd.z_bed, :, :, 1)
    shakti_sim.state.zs     .= interior(yelmo.tpo.z_srf, :, :, 1)
    shakti_sim.state.H      .= interior(yelmo.tpo.H_ice, :, :, 1)
    shakti_sim.state.A_visc .= perYear2perSecond.(mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1])
    shakti_sim.state.ub_x   .= perYear2perSecond.(interior(yelmo.dyn.ux_b, :, :, 1))
    shakti_sim.state.ub_y   .= perYear2perSecond.(interior(yelmo.dyn.uy_b, :, :, 1))
    shakti_sim.state.taub_x .= interior(yelmo.dyn.taub_acx, :, :, 1)
    shakti_sim.state.taub_y .= interior(yelmo.dyn.taub_acy, :, :, 1)
    shakti_sim.state.G      .= interior(yelmo.bnd.Q_geo, :, :, 1) .* 1e-3
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
    Yelmo_to_FastHydrology_Shakti!(shakti_sim, yelmo)
    # dt is in years (Yelmo's unit); shakti_sim.dt is in seconds -- sub-cycle enough Shakti
    # steps to cover one Yelmo timestep.
    dt_seconds = dt * FastHydrology.SECONDS_PER_YEAR
    n_sub = min(SHAKTI_MAX_SUBSTEPS_SMOKE_TEST, max(1, round(Int, dt_seconds / shakti_sim.dt)))
    for _ in 1:n_sub
        FastHydrology.step!(sim)
    end
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

"""Build a fresh, initialized Yelmo model (same recipe every time, so runs start from the same t=0 state)."""
function setup_yelmo()
    p, yelmo = build_yelmo()
    yelmo.bnd.H_sed .= 100.0
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

    plot_title_list    = ["No coupling", "K24 coupling", "HAB coupling", "Shakti coupling"]
    make_coupling_list = [
        yelmo -> NoCoupling(),
        yelmo -> CoupledHydrology(build_hydrology_sim_K24(yelmo)),
        yelmo -> CoupledHydrology(build_hydrology_sim_HAB(yelmo)),
        yelmo -> CoupledHydrology(build_hydrology_sim_Shakti(yelmo)),
    ]

    for (idx, make_coupling) in enumerate(make_coupling_list)

        # Fresh Yelmo model per branch, so all four start from the same t=0 state --
        # otherwise each branch would continue from wherever the previous one left off.
        yelmo    = setup_yelmo()
        coupling = make_coupling(yelmo)

        # Write to file
        # yelmo_output_groups = [:tpo, :dyn, :thrm, :mat, :bnd]
        # yelmo_out = init_output(yelmo, joinpath(yelmo.rundir, "yelmo.nc"); selection=OutputSelection(groups=yelmo_output_groups))
        # write_output!(yelmo_out, yelmo)

        # Run the total time simulation
        @time run!(coupling, yelmo; dt=2.0, time_end=4.0)

        # close(yelmo_out)

        # Visualize fields -- headless-safe: saved to PLOT_DIR rather than displayed, since a
        # cluster run has no display to pop a window on.
        slug = lowercase(replace(plot_title_list[idx], " " => "_"))
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
