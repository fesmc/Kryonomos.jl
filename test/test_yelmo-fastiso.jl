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

# Define parameters and write parameter file
p = YelmoParameters("Greenland")

# Initialize Yelmo
ylmo = YelmoMirror(p, 0.0; rundir="run01", overwrite=true);

# Initialize FastIsostasy #################
W, n, T = 3f6, 7, Float32
Nx, Ny = ylmo.g.Nx, ylmo.g.Ny
domain = RegionalDomain(3500e3, 3500e3, Nx, Ny) 
c = PhysicalConstants(rho_ice = ylmo.p.phys.rho_ice)

bcs = BoundaryConditions(domain, ice_thickness = ExternallyUpdatedIceThickness)
sealevel = RegionalSeaLevel(
    surface = LaterallyVariableSeaSurface(),
    load = InteractiveSealevelLoad(),
    bsl = PiecewiseConstantBSL(),
)

G = 0.50605f11              # Shear modulus (Pa)
nu = 0.28f0                 # Poisson ratio
E = G * 2 * (1 + nu)        # Young modulus
lb = c.r_equator .- [6301f3, 5951f3, 5701f3]    # 3 layer boundaries
solidearth = SolidEarth(
    domain,
    layer_boundaries = T.(lb),
    layer_viscosities = [1f21, 1f21, 2f21],
    litho_youngmodulus = E,
    litho_poissonratio = nu,
    rho_litho = 3100f0,
    rho_uppermantle = 3500f0,
)

###########################################

# Populate boundary fields
ylmo.bnd.H_sed .= 100.0;

# Initialize Yelmo state
init_state!(ylmo, 0.0, "robin-cold");

time_init, time_end, dt = 0.0, 5.0, 1.0;

for t in time_init:dt:time_end

    # Update isostasy
    ylmo.bnd.z_bed .+= 100.0

    # Advance by dt
    time_step!(ylmo,t-ylmo.time);

    write_output!(out, ylmo)
end

# Close output file
close(out)

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
