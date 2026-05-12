## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using Revise
using CairoMakie
using Oceananigans.Fields
using Oceananigans.Grids
using Yelmo
using FastIsostasy

# Initialize Yelmo #################

p = YelmoParameters("Greenland")
ylmo = YelmoMirror(p, 0.0; rundir="run01", overwrite=true);

# Populate boundary fields
ylmo.bnd.H_sed .= 100.0;

# Initialize Yelmo state
init_state!(ylmo, 0.0; thrm_method="robin-cold");

# Initialize FastIsostasy ##########
# Identical-grid coupling: FastIsostasy uses the same (Nx, Ny, dx, dy)
# as Yelmo so push/pull is a plain broadcast copy. To run FastIsostasy
# on a different grid, swap `push_ice!` / `pull_bedrock!` below for
# interpolating versions.

T = Float32
time_init, time_end, dt = 0.0, 5.0, 1.0
t_out = collect(T(time_init):T(dt):T(time_end))

Nx, Ny = ylmo.g.Nx, ylmo.g.Ny
Wx, Wy = T(ylmo.g.Δxᶜᵃᵃ * Nx), T(ylmo.g.Δyᵃᶜᵃ * Ny)
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
# grid, swap these for interpolating versions).
push_ice!(sim, ylmo) = (sim.now.H_ice .= ylmo.tpo.H_ice; nothing)
pull_bedrock!(ylmo, sim) = begin
    ylmo.bnd.z_bed .= ylmo.bnd.z_bed_ref .+ sim.now.u .+ sim.now.ue
    ylmo.bnd.z_sl  .= sim.now.z_ss
    return nothing
end

# Initialize Yelmo output file ####
yelmo_out = init_output(ylmo, joinpath(ylmo.rundir, "yelmo.nc"),
    selection = OutputSelection(
        groups = [:tpo, :dyn, :thrm, :mat, :bnd],
    )
)
write_output!(yelmo_out, ylmo)

# Time loop (Yelmo only — coupling lands in commit 4)
for t in T(time_init):T(dt):T(time_end)
    Yelmo.step!(ylmo, t - ylmo.time)
    write_output!(yelmo_out, ylmo)
end

close(yelmo_out)

# Plot some data
heatmap(ylmo.dyn.uxy_s, colorscale = log10)

# Standalone FastIsostasy smoke test — step the integrator forward by
# `dt` repeatedly with H_ice held at zero. Confirms FastIsostasy
# initializes and integrates cleanly before we couple it in commit 4.
println("FastIsostasy smoke test:")
println("  initial extrema(u)  = ", extrema(sim.now.u))
println("  initial extrema(ue) = ", extrema(sim.now.ue))
for t in T(time_init):T(dt):T(time_end - dt)
    FastIsostasy.step!(integrator, T(dt), true)
    println("  t=$(t + dt)  extrema(u)=$(extrema(sim.now.u))  extrema(ue)=$(extrema(sim.now.ue))")
end
