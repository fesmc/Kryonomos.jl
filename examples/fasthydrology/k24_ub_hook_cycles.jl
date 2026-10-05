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

    fortran_k24 = get(ENV, "YELMO_YHYD_K24", "0") == "1"   # no Julia hydrology: Fortran Yelmo runs K24 (and its own hook)
    yelmo    = get(ENV, "WARM", "0") == "1" ? setup_warm() : setup_yelmo(; external_neff = !fortran_k24)
    coupling = fortran_k24 ? NoCoupling() : CoupledHydrology(build_hydrology_sim_K24(yelmo))
    install_ub_hook!(coupling, yelmo)
    @info "hook installed" installed = (yelmo.hooks.neff_from_ub !== nothing) fortran_k24

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
        Nprev = copy(interior(yelmo.dyn.N_eff, :, :, 1)); NEFF_IN[] = nothing
        step!(coupling, yelmo, k * dt, dt)
        if NEFF_IN[] !== nothing   # N the step started with (full K24 update) vs the previous step's final N
            gm = interior(yelmo.tpo.f_grnd, :, :, 1) .> 0.5
            r = abs.(NEFF_IN[] .- Nprev)[gm] ./ max.(abs.(Nprev[gm]), 1e3)
            @info "  step-start N vs previous final N: cells off by >10%: $(count(>(0.1), r)) of $(length(r))"
        end
        hasproperty(yelmo.dyn, :scratch) && push!(iters, yelmo.dyn.scratch.ssa_iter_now[])
        save(k + 1)
        if coupling isa CoupledHydrology && coupling.sim.model isa KazmierczakHydroModel   # is N the one of the step's final u_b? (routing held)
            sim = coupling.sim
            Nout = copy(interior(yelmo.dyn.N_eff, :, :, 1))
            Nchk = copy(interior(FastHydrology.N_from_ub!(sim.model, sim.grid, sim.state,
                       perYear2perSecond.(interior(yelmo.dyn.uxy_b, :, :, 1))), :, :, 1))
            gm = interior(yelmo.tpo.f_grnd, :, :, 1) .> 0.5
            rel = abs.(Nout .- Nchk)[gm] ./ max.(abs.(Nchk[gm]), 1e3)
            @info "  N(final u_b) consistency: cells off by >10%: $(count(>(0.1), rel)) of $(length(rel)), median rel $(round(median(rel); sigdigits = 2))"
            if get(ENV, "K24_DECOMPOSE", "0") == "1"   # what changes in the routing between consecutive full updates?
                m = sim.model; L = FastHydrology.KAZMIERCZAK_DEFAULT_L_W
                old = (fixed = copy(interior(m.mdot_fixed, :, :, 1)), Qb = copy(interior(m.Q_b, :, :, 1)),
                       Qd = copy(interior(m.Q_diss, :, :, 1)), tot = copy(interior(m.mdot_total, :, :, 1)),
                       q = copy(interior(m.q, :, :, 1)), ub = copy(interior(m.abs_v_b, :, :, 1)))
                Yelmo_to_FastHydrology_K24!(sim, yelmo); FastHydrology.run!(sim)
                Nfull = copy(interior(sim.state.N, :, :, 1))
                bad = gm .& (abs.(Nfull .- Nchk) ./ max.(abs.(Nchk), 1e3) .> 0.1)
                med(x) = isempty(x) ? NaN : round(median(x); sigdigits = 3)
                nw = (fixed = interior(m.mdot_fixed, :, :, 1), Qb = interior(m.Q_b, :, :, 1), Qd = interior(m.Q_diss, :, :, 1),
                      tot = interior(m.mdot_total, :, :, 1), q = interior(m.q, :, :, 1))
                @info "  decompose ($(count(bad)) cells where full-update N differs >10% from held-routing N), medians over them:" *
                      " mdot_total old/new $(med(old.tot[bad])) / $(med(nw.tot[bad]))" *
                      "  |dfixed|/|tot| $(med(abs.(nw.fixed[bad] .- old.fixed[bad]) ./ max.(abs.(nw.tot[bad]), 1e-14)))" *
                      "  |dQb/L|/|tot| $(med(abs.(nw.Qb[bad] .- old.Qb[bad]) ./ L ./ max.(abs.(nw.tot[bad]), 1e-14)))" *
                      "  |dQdiss/L|/|tot| $(med(abs.(nw.Qd[bad] .- old.Qd[bad]) ./ L ./ max.(abs.(nw.tot[bad]), 1e-14)))" *
                      "  |dq|/q $(med(abs.(nw.q[bad] .- old.q[bad]) ./ max.(nw.q[bad], 1e-14)))"
            end
        end
        @info "step $k  t=$(round(yelmo.time; digits = 3)) yr  hook calls $(NEFF_CALLS[])  Picard iters $(isempty(iters) ? "-" : iters[end])  $(round(time() - t0; digits = 1)) s"
    end
    close(ds)
    isempty(iters) || @info "mean Picard iterations per step" mean(iters)
end

main_cycles(ARGS[1], parse(Float64, get(ARGS, 2, "0.1")), parse(Float64, get(ARGS, 3, "2.0")))
