# Short K24-coupled Greenland run on the Julia Yelmo backend, writing N_eff and f_grnd at every
# step to <outdir>/yelmo.nc so that ~/build_jobs/cyc.jl can count the cells whose N alternates by
# more than 10% over the last 7 outputs. Run it twice, K24_UB_HOOK=0 and K24_UB_HOOK=1:
#
#   julia --project=<examples env> k24_ub_hook_cycles.jl <outdir> [dt_yr=0.1] [time_end_yr=2.0]
#
# Reference (Fortran, GRL-16KM, 2 yr in 0.1 yr steps): 3215 cycling cells without the hook, 95-562 with it.
ENV["KRYO_NO_MAIN"] = "1"
include(joinpath(@__DIR__, "Greenland_yelmo-fasthydrology.jl"))

# WARM=1: start from the ISMIP7 15 kyr GRL-16KM spin-up restart (optimised friction, spun-up
# temperature) instead of the cold InitMIP start. Parameters from the ISMIP7 namelist as the Julia
# structs read it; N_eff is set externally (yneff.method = -1) and cb_ref is kept from the restart.
function setup_warm()
    p = YelmoParameters(YELMO_NML_ISMIP7, "Greenland")
    p = _override_field(p, :yneff, _override_field(p.yneff, :method, -1))
    mkpath(RUN_DIR)
    yelmo = YelmoModel(ISMIP7_RESTART, 0.0; p, boundaries = :bounded, rundir = RUN_DIR, strict = false)
    @info "warm start" restart = ISMIP7_RESTART ytill_method = p.ytill.method solver = p.ydyn.solver
    Yelmo.update_diagnostics!(yelmo)
    Yelmo.YelmoModelDyn.dyn_step!(yelmo, 0.0)
    return yelmo
end

function main_cycles(outdir, dt, time_end)
    mkpath(outdir)
        @info "K24 on Yelmo ($BACKEND)" hook = K24_UB_HOOK sliding = K24_SLIDING dt time_end

    yelmo    = get(ENV, "WARM", "0") == "1" ? setup_warm() : setup_yelmo(; external_neff = true)
    coupling = CoupledHydrology(build_hydrology_sim_K24(yelmo))
    install_ub_hook!(coupling, yelmo)
    @info "hook installed" installed = (yelmo.hooks.neff_from_ub !== nothing)

    nx, ny = yelmo.g.Nx, yelmo.g.Ny
    nsteps = round(Int, time_end / dt)
    ds = NCDataset(joinpath(outdir, "yelmo.nc"), "c")
    defDim(ds, "x", nx); defDim(ds, "y", ny); defDim(ds, "time", nsteps + 1)
    vt = defVar(ds, "time", Float64, ("time",))
    vN = defVar(ds, "N_eff", Float64, ("x", "y", "time"))
    vg = defVar(ds, "f_grnd", Float64, ("x", "y", "time"))
    vu = defVar(ds, "uxy_b", Float64, ("x", "y", "time"))
    iters = Int[]
    function save(k)
        vt[k]        = yelmo.time
        vN[:, :, k]  = interior(yelmo.dyn.N_eff, :, :, 1)
        vg[:, :, k]  = interior(yelmo.tpo.f_grnd, :, :, 1)
        vu[:, :, k]  = interior(yelmo.dyn.uxy_b, :, :, 1)
    end
    save(1)
    for k in 1:nsteps
        t0 = time()
        step!(coupling, yelmo, k * dt, dt)
        hasproperty(yelmo.dyn, :scratch) && push!(iters, yelmo.dyn.scratch.ssa_iter_now[])
        save(k + 1)
        if coupling.sim.model isa KazmierczakHydroModel   # is N the one of the step's final u_b? (routing held)
            sim = coupling.sim
            Nout = copy(interior(yelmo.dyn.N_eff, :, :, 1))
            Nchk = copy(interior(FastHydrology.N_from_ub!(sim.model, sim.grid, sim.state,
                       perYear2perSecond.(interior(yelmo.dyn.uxy_b, :, :, 1))), :, :, 1))
            gm = interior(yelmo.tpo.f_grnd, :, :, 1) .> 0.5
            rel = abs.(Nout .- Nchk)[gm] ./ max.(abs.(Nchk[gm]), 1e3)
            @info "  N(final u_b) consistency: cells off by >10%: $(count(>(0.1), rel)) of $(length(rel)), median rel $(round(median(rel); sigdigits = 2))"
        end
        @info "step $k  t=$(round(yelmo.time; digits = 3)) yr  hook calls $(NEFF_CALLS[])  Picard iters $(isempty(iters) ? "-" : iters[end])  $(round(time() - t0; digits = 1)) s"
    end
    close(ds)
    isempty(iters) || @info "mean Picard iterations per step" mean(iters)
end

main_cycles(ARGS[1], parse(Float64, get(ARGS, 2, "0.1")), parse(Float64, get(ARGS, 3, "2.0")))
