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
using FastIsostasy
using IceSheetBenchmarks
using Statistics
using FastHydrology

# ── Constants ─────────────────────────────────────────────────────────────────

const YELMO_PROJECT                = "/Users/taange001/Documents/Coding/Yelmo.jl"
const ICE_SHEET_BENCHMARKS_PROJECT = joinpath(YELMO_PROJECT, "benchmarks/IceSheetBenchmarks")
const RUN_DIR                      = @__DIR__
const DATA_DIR                     = joinpath(YELMO_PROJECT, "benchmarks/initmip-grl/data/GRL-16KM")

const YELMO_NML          = "/Users/taange001/Documents/Coding/Kryonomos.jl/examples/run01/Greenland.nml"
const YELMO_REGIONS_FILE = joinpath(DATA_DIR, "GRL-16KM_REGIONS.nc")
const YELMO_BASINS_FILE  = joinpath(DATA_DIR, "GRL-16KM_BASINS-nasa.nc")
const YELMO_TOPO_FILE    = joinpath(DATA_DIR, "GRL-16KM_TOPO-M17-v5.nc")

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

"""Fully coupled: Yelmo ↔ FastHydrology exchange each timestep."""
struct CoupledHydrology <: HydrologyCoupling
    sim::SteadyStateSimulation
end

"""No hydrology: Yelmo advances alone."""
struct NoCoupling <: HydrologyCoupling end

# ── Coupling functions ────────────────────────────────────────────────────────

"""Copy Yelmo fields → FastHydrology state/model."""
function Yelmo_to_FastHydrology!(sim, yelmo)
    sim.state.mask    .= _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd) .> 0
    sim.state.h       .= interior(yelmo.tpo.H_ice, :, :, 1)
    sim.state.b       .= interior(yelmo.bnd.z_bed, :, :, 1)
    sim.model.abs_v_b .= interior(yelmo.dyn.uxy_b, :, :, 1)
    sim.model.A_visc  .= mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1]
    sim.model.mdot    .= perYear2perSecond.(interior(yelmo.thrm.bmb_grnd, :, :, 1))
    sim.model.kappa   .= zeros(T, Nx, Ny)
end

"""Copy FastHydrology outputs → Yelmo boundary fields."""
function FastHydrology_to_Yelmo!(yelmo, sim)
    yelmo.dyn.N_eff .= sim.state.N   # effective pressure
    yelmo.thrm.H_w  .= sim.state.W   # subglacial water thickness
end

# ── Time stepping (dispatched on coupling type) ───────────────────────────────

"""Single timestep with active hydrology coupling."""
function step!(coupling::CoupledHydrology, yelmo, t, dt)
    Yelmo_to_FastHydrology!(coupling.sim, yelmo)
    FastHydrology.run!(coupling.sim)
    FastHydrology_to_Yelmo!(yelmo, coupling.sim)
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
        # println("t = $t yr : $(extrema(yelmo.dyn.uxy_s))")
    end
end

function main()

    # ── Yelmo setup ───────────────────────────────────────────────────────────────

    p, yelmo = build_yelmo()
    yelmo.bnd.H_sed .= 100.0
    init_state!(yelmo, 0.0; thrm_method="robin-cold")

    # ── FastHydrology setup ───────────────────────────────────────────────────────

    T = Float32

    Nx, Ny = yelmo.g.Nx, yelmo.g.Ny
    xlims  = (0, T(yelmo.g.Δxᶜᵃᵃ * Nx))
    ylims  = (0, T(yelmo.g.Δyᵃᶜᵃ * Ny))

    # Initial fields from Yelmo
    mask    = _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd) .> 0
    h       = interior(yelmo.tpo.H_ice,  :, :, 1)   # ice thickness
    b       = interior(yelmo.bnd.z_bed,  :, :, 1)   # bedrock elevation
    abs_v_b = interior(yelmo.dyn.uxy_b,  :, :, 1)   # basal speed
    A_visc  = mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1]  # depth-averaged rate factor
    mdot    = perYear2perSecond.(yelmo.thrm.bmb_grnd)          # basal melt rate [kg m⁻² s⁻¹]
    kappa   = zeros(T, Nx, Ny)                                 # bed hardness (0: hard, 1: soft)

    longcoupwater = 0.0   # smoothing of geometric-potential gradients
    fill_iters    = 10    # iterations to fill local minima in potential field

    grid  = OGRectHydroGrid(yelmo.g)
    model = KazmierczakHydroModel(grid, kappa, abs_v_b, A_visc, mdot; longcoupwater=longcoupwater, fill_iters=fill_iters)
    state = HydroState(grid, mask, h, b)
    sim   = SteadyStateSimulation(model, grid, state)

    # ── Run ───────────────────────────────────────────────────────────────────────

    plot_title_list = ["No coupling", "Hydrology coupling"]

    for (idx, coupling) in enumerate([NoCoupling(), CoupledHydrology(sim)])

        # Write to file
        # yelmo_output_groups = [:tpo, :dyn, :thrm, :mat, :bnd]
        # yelmo_out = init_output(yelmo, joinpath(yelmo.rundir, "yelmo.nc"); selection=OutputSelection(groups=yelmo_output_groups))
        # write_output!(yelmo_out, yelmo)

        # Run the total time simulation
        @time run!(coupling, yelmo; dt=2.0, time_end=50.0)

        # close(yelmo_out)

        # Visualize fields
        # visualize_field(model.q; plot_title="Water flux")
        visualize_field(state.W; plot_title="Water thickness: " * plot_title_list[idx])
        visualize_field(yelmo.dyn.uxy_s; plot_title = "Surface horizontal velocity magnitude: " * plot_title_list[idx])

    end

end

main()
