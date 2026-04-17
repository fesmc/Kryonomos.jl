## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using BenchmarkTools
using Revise
using CairoMakie
using Oceananigans.Fields
using Oceananigans.Grids
using Yelmo
using FastIsostasy

# Initialize Yelmo #################

# Define parameters and write parameter file
p = YelmoParameters("Greenland")

# Initialize Yelmo
ylmo = YelmoMirror(p, 0.0; rundir="run01", overwrite=true);

# Initialize FastIsostasy #################

W, n, T = 3f6, 7, Float32
Nx, Ny = ylmo.g.Nx, ylmo.g.Ny
domain = RegionalDomain(T(3500e3), T(3500e3), Nx, Ny) 
c = PhysicalConstants(rho_ice = ylmo.p.phys.rho_ice)

bcs = BoundaryConditions(domain, ice_thickness = ExternallyUpdatedIceThickness())
sealevel = RegionalSeaLevel(
    surface = LaterallyVariableSeaSurface(),
    load = NoSealevelLoad(),
    bsl = PiecewiseConstantBSL(),
)

solidearth = SolidEarth(
    domain,
)
solidearth.litho_rigidity .= T(1f25)
solidearth.effective_viscosity .= T(1f21)

sim1 = Simulation(domain, bcs, sealevel, solidearth, (0, 1f3))
integrator = init_integrator(sim1::Simulation)

###########################################

# Populate boundary fields
ylmo.bnd.H_sed .= 100.0;

# Initialize Yelmo state
init_state!(ylmo, 0.0, "robin-cold");

# Initialize output files
yelmo_out = init_output(ylmo, joinpath(ylmo.rundir,"yelmo.nc"),
    selection = OutputSelection(
        groups  = [:tpo,:dyn,:thrm,:mat,:bnd],
    )
)
time_init, time_end, dt = 0.0, 5.0, 1.0;

for t in time_init:dt:time_end

    # Update isostasy
    isos.now.H_ice .= 
    step!()
    ylmo.bnd.z_bed .= ylmo.bnd.z_ref .+ isos.now.u .+ isos.now.ue
    ylmo.bnd.z_sl .= isos.now.z_ss

    # Advance by dt
    time_step!(ylmo,t-ylmo.time);

    write_output!(yelmo_out, ylmo)
end

# Close output file
close(yelmo_out)

# Plot some data
heatmap(ylmo.dyn.uxy_s,colorscale=log10)



# Initialize output file
out = init_output(ylmo, joinpath(ylmo.rundir,"yelmo.nc"),
    selection = OutputSelection(
        groups  = [:tpo,:dyn,:thrm,:mat,:bnd],
    )
)

# Write initial state
write_output!(out, ylmo)
#close(out)
