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

const YELMO_PROJECT            = "/Users/taange001/Documents/Coding/Yelmo.jl"
const ICE_SHEET_BENCHMARKS_PROJECT = joinpath(YELMO_PROJECT, "benchmarks/IceSheetBenchmarks")
const RUN_DIR                  = @__DIR__
const DATA_DIR                 = joinpath(YELMO_PROJECT, "benchmarks/initmip-grl/data/GRL-16KM")

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
    masks = _override_field(masks,         :basins_path,   YELMO_BASINS_FILE)

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

p, yelmo = build_yelmo()

yelmo.bnd.H_sed .= 100.0
init_state!(yelmo, 0.0; thrm_method="robin-cold")

# ── FastHydrology setup ───────────────────────────────────────────────────────

const T = Float32

time_init, time_end, dt = 0.0, 1.0, 1.0   # years

Nx, Ny = yelmo.g.Nx, yelmo.g.Ny
xlims  = (0, T(yelmo.g.Δxᶜᵃᵃ * Nx))
ylims  = (0, T(yelmo.g.Δyᵃᶜᵃ * Ny))

# Initial fields from Yelmo
mask      = _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd) .> 0
h         = interior(yelmo.tpo.H_ice,    :, :, 1)   # ice thickness
b         = interior(yelmo.bnd.z_bed,    :, :, 1)   # bedrock elevation
abs_v_b   = interior(yelmo.dyn.uxy_b,   :, :, 1)   # basal speed
A_visc    = mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1]   # depth-averaged rate factor
mdot      = perYear2perSecond.(yelmo.thrm.bmb_grnd)          # basal melt rate [kg m⁻² s⁻¹]
kappa     = zeros(T, Nx, Ny)                                  # bed hardness (0: hard, 1: soft)

const LONGCOUPWATER = 0.0   # smoothing of geometric-potential gradients
const FILL_ITERS    = 10    # iterations to fill local minima in potential field

grid  = OGRectHydroGrid(yelmo.g)
model = KazmierczakHydroModel(grid, kappa, abs_v_b, A_visc, mdot;
                               longcoupwater=LONGCOUPWATER, fill_iters=FILL_ITERS)
state = HydroState(grid, mask, h, b)
sim   = SteadyStateSimulation(model, grid, state)

# ── Coupling functions ────────────────────────────────────────────────────────

"""Copy Yelmo fields → FastHydrology state/model."""
function yelmo_to_FastHydrology!(sim, yelmo)
    sim.state.mask  .= _compute_interior(yelmo.tpo.f_ice * yelmo.tpo.f_grnd) .> 0
    sim.state.h     .= interior(yelmo.tpo.H_ice,  :, :, 1)
    sim.state.b     .= interior(yelmo.bnd.z_bed,  :, :, 1)

    sim.model.abs_v_b .= interior(yelmo.dyn.uxy_b, :, :, 1)
    sim.model.A_visc  .= mean(interior(yelmo.mat.ATT), dims=3)[:, :, 1]
    sim.model.mdot    .= perYear2perSecond.(interior(yelmo.thrm.bmb_grnd, :, :, 1))
    sim.model.kappa   .= zeros(T, Nx, Ny)
end

"""Copy FastHydrology outputs → Yelmo boundary fields."""
function FastHydrology_to_yelmo!(yelmo, sim)
    yelmo.dyn.N_eff .= sim.state.N   # effective pressure
    yelmo.thrm.H_w  .= sim.state.W   # subglacial water thickness
end

# ── Output ────────────────────────────────────────────────────────────────────

const YELMO_OUTPUT_GROUPS = [:tpo, :dyn, :thrm, :mat, :bnd]

yelmo_out = init_output(yelmo, joinpath(yelmo.rundir, "yelmo.nc"); selection=OutputSelection(groups=YELMO_OUTPUT_GROUPS))
write_output!(yelmo_out, yelmo)

# ── Coupled time loop ─────────────────────────────────────────────────────────

function run!()
    for t in dt:dt:time_end
        # yelmo_to_FastHydrology!(sim, yelmo)
        # FastHydrology.run!(sim)
        # FastHydrology_to_yelmo!(yelmo, sim)
        Yelmo.step!(yelmo, dt)
        write_output!(yelmo_out, yelmo)
        println("t = $t yr")
    end
end

@time run!()
close(yelmo_out)

# ── Diagnostics ───────────────────────────────────────────────────────────────

fig = Figure()
ax  = Axis(fig[1, 1])
hm  = heatmap!(ax, interior(yelmo.thrm.bmb_grnd, :, :, 1))
Colorbar(fig[1, 2], hm)
display(fig)