## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

using Revise
using CairoMakie
using Oceananigans.Fields
using Oceananigans.Grids
using Yelmo

# Initialize Yelmo #################

p = YelmoParameters("Greenland")
ylmo = YelmoMirror(p, 0.0; rundir="run01", overwrite=true);

# Populate boundary fields
ylmo.bnd.H_sed .= 100.0;

# Initialize Yelmo state
init_state!(ylmo, 0.0; thrm_method="robin-cold");

# Initialize output file
yelmo_out = init_output(ylmo, joinpath(ylmo.rundir, "yelmo.nc"),
    selection = OutputSelection(
        groups = [:tpo, :dyn, :thrm, :mat, :bnd],
    )
)
write_output!(yelmo_out, ylmo)

# Time loop (Yelmo only — no isostasy yet)
time_init, time_end, dt = 0.0, 5.0, 1.0;

for t in time_init:dt:time_end
    step!(ylmo, t - ylmo.time)
    write_output!(yelmo_out, ylmo)
end

close(yelmo_out)

# Plot some data
heatmap(ylmo.dyn.uxy_s, colorscale = log10)
