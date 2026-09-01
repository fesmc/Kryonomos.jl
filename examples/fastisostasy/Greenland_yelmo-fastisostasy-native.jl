## Preamble #############################################
# Run from examples/: that is where Project.toml/Manifest.toml live, and where the
# Fortran side resolves its relative input/ paths and ice_data/. This script lives
# one level down, so cd to the parent rather than to @__DIR__.
cd(dirname(@__DIR__))
import Pkg; Pkg.activate(".")
#########################################################

# Coupled Yelmo + FastIsostasy example using the *native* (pure-Julia)
# `YelmoModel` instead of the Fortran-backed `YelmoMirror`. The model is
# built from scratch on the GRL-16KM grid (read from a NetCDF grid file),
# with its vertical axes derived from the Yelmo parameters, and its
# initial + reference topography loaded from a TOPO file via the
# `init_topo_load!` routine.

using Revise
using CairoMakie
using Oceananigans.Fields
using Oceananigans.Grids
using Yelmo
using FastIsostasy

# Initialize Yelmo (native YelmoModel) #################

# Default parameters throughout; only the initial-topography file is
# pointed at the M17-v5 variant. `domain`/`grid_name` default to
# "Greenland"/"GRL-16KM", so the {domain}/{grid_name} template resolves
# to the GRL-16KM ice_data files.
# Note: `yelmo_init_topo_params` is exported by both the mirror (`YelmoPar`)
# and native (`YelmoModelPar`) parameter modules, so qualify the native one.
p = YelmoModelParameters("Greenland";
    # Use the "fixed" dynamics + thermodynamics solvers. The native ports
    # of the `diva` velocity solver and the `robin` thermodynamics solver
    # currently blow up (H_ice / thermal state -> NaN) after one 1-yr step
    # on the real GRL-16KM topography, so we hold ice velocity and the
    # thermal state fixed. The ice geometry then stays put while the
    # bedrock responds to the (static) Greenland ice load via FastIsostasy
    # — a clean demonstration of the Yelmo<->FastIsostasy coupling.
    ydyn   = Yelmo.YelmoModelPar.ydyn_params(solver = "fixed"),
    ytherm = Yelmo.YelmoModelPar.ytherm_params(method = "fixed"),
    yelmo_init_topo = Yelmo.YelmoModelPar.yelmo_init_topo_params(
        init_topo_path = "ice_data/{domain}/{grid_name}/{grid_name}_TOPO-M17-v5.nc",
    ),
)

# Build the model from scratch on the grid defined in the REGIONS file.
# Only the (xc, yc) axes are read here; vertical axes come from `p`, and
# all state fields start at their default (zero) allocation.
gridfile = "ice_data/Greenland/GRL-16KM/GRL-16KM_REGIONS.nc"
y = YelmoModel(gridfile, p; rundir="run01", alias="ymodel")

# Populate initial topography (H_ice, z_bed, z_bed_sd, z_srf) from the
# TOPO-M17-v5 file, then set the sediment field as in the mirror example.
init_topo_load!(y)
interior(y.bnd.H_sed) .= 100.0

# Reference bedrock for isostasy: snapshot the loaded (initial) z_bed.
# FastIsostasy displacements (u + ue) are applied relative to this.
z_bed_ref = copy(interior(y.bnd.z_bed)[:, :, 1])

# Initialize Yelmo state
init_state!(y, 0.0; thrm_method="robin-cold")

# Initialize FastIsostasy ##########
# Identical-grid coupling: FastIsostasy uses the same (Nx, Ny, dx, dy)
# as Yelmo so push/pull is a plain broadcast copy. To run FastIsostasy
# on a different grid, swap `push_ice!` / `pull_bedrock!` below for
# interpolating versions.

T = Float32
time_init, time_end, dt = 0.0, 5.0, 1.0
t_out = collect(T(time_init):T(dt):T(time_end))

Nx, Ny = y.g.Nx, y.g.Ny
# RegionalDomain treats (Wx, Wy) as half-widths (it builds a grid on
# [-Wx, Wx] with dx = 2·Wx/Nx). To match Yelmo's grid spacing exactly,
# pass half the total extent so FastIsostasy's dx equals Yelmo's.
Wx, Wy = T(y.g.Δxᶜᵃᵃ * Nx / 2), T(y.g.Δyᵃᶜᵃ * Ny / 2)
domain = RegionalDomain(Wx, Wy, Nx, Ny)
bcs = BoundaryConditions(domain, ice_thickness = ExternallyUpdatedIceThickness())
sealevel = RegionalSeaLevel(
    surface = LaterallyVariableSeaSurface(),
    load    = NoSealevelLoad(),
    bsl     = PiecewiseConstantBSL(),
)
solidearth = SolidEarth(domain)
sim = Simulation(domain, bcs, sealevel, solidearth,
    (T(time_init), T(time_end));
    nout = NativeOutput(t = t_out),
)
integrator = init_integrator(sim)

# Coupling helpers (identical-grid case — for a different FastIsostasy
# grid, swap these for interpolating versions). The native YelmoModel
# stores 2D fields as Oceananigans Center fields, so we couple through
# their interior (dropping the singleton vertical dimension).
push_ice!(sim, y) = (sim.now.H_ice .= @view interior(y.tpo.H_ice)[:, :, 1]; nothing)
pull_bedrock!(y, sim, z_bed_ref) = begin
    @view(interior(y.bnd.z_bed)[:, :, 1]) .= z_bed_ref .+ sim.now.u .+ sim.now.ue
    @view(interior(y.bnd.z_sl)[:, :, 1])  .= sim.now.z_ss
    return nothing
end

# Initialize Yelmo output file ####
yelmo_out = init_output(y, joinpath(y.rundir, "yelmo.nc"),
    selection = OutputSelection(
        groups = [:tpo, :dyn, :thrm, :mat, :bnd],
    )
)
write_output!(yelmo_out, y)

# Coupled time loop
for t in dt:dt:time_end
    push_ice!(sim, y)
    FastIsostasy.step!(integrator, T(dt), true)
    pull_bedrock!(y, sim, z_bed_ref)
    Yelmo.step!(y, dt)
    write_output!(yelmo_out, y)
    println("t=$t  extrema(u)=$(extrema(sim.now.u))  extrema(z_bed)=$(extrema(interior(y.bnd.z_bed)))")
end

close(yelmo_out)

# Plot some data
heatmap(interior(y.dyn.uxy_s)[:, :, 1], colorscale = log10)
