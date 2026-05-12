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
# Standalone setup (no coupling to Yelmo yet — grid is hardcoded).
# Commit 3 will derive the domain from Yelmo's grid.

T = Float32
time_init, time_end, dt = 0.0, 5.0, 1.0
t_out = collect(T(time_init):T(dt):T(time_end))

domain = RegionalDomain(T(3500e3), T(3500e3), 64, 64)
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
