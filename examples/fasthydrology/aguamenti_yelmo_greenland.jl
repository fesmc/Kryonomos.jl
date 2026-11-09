# Aguamenti v1.0 (FastHydrology.jl) coupled to Yelmo (Julia backend) on Greenland 16 km, from the
# ISMIP7 15 kyr spin-up restart (WARM=1, default) or the cold InitMIP start (WARM=0).
#
#   julia --project=<env> aguamenti_yelmo_greenland.jl <outdir> [dt_yr=0.1] [time_end_yr=2.0]
#
# Env: HYDRO_MODEL = aguamenti (default) | k24 (the same run with K24, for comparison);
#      K24_UB_HOOK = 1 (default: N re-evaluated in every DIVA iteration from that iteration's u_b) | 0;
#      AGU_KAPPA, AGU_COLD_T, AGU_SPINUP (see Greenland_yelmo-fasthydrology.jl).
# Writes N_eff, f_grnd, uxy_b (and for Aguamenti regime, chi, N_t, dry) at every step to
# <outdir>/yelmo.nc; ~/build_jobs/cyc.jl counts the cells whose N alternates by more than 10 %.
ENV["KRYO_NO_MAIN"] = "1"
include(joinpath(@__DIR__, "Greenland_yelmo-fasthydrology.jl"))

const HYDRO_MODEL = get(ENV, "HYDRO_MODEL", "aguamenti")

# The ISMIP7 namelist is written for Fortran Yelmo (YelmoMirror). For the Julia backend its `slide_T`
# (&ydyn, no temperature scaling of sliding) is `frz_scale`, and keys Yelmo.jl does not know
# (`use_strain_sia` in &ytherm) are dropped with a warning (read_nml(...; strict = false)).
# This driver runs the Julia backend only (FASTHYDRO_BACKEND=yelmo).
function ismip7_parameters_julia()
    txt = replace(read(YELMO_NML_ISMIP7, String), r"(?m)^(\s*)slide_T(\s*)=" => s"\1frz_scale\2=")
    mkpath(RUN_DIR)
    path = joinpath(RUN_DIR, "Greenland_ismip7_julia.nml")
    write(path, txt)
    P = Yelmo.YelmoPar
    p0 = P.read_nml(path; strict = false)
    # Options YelmoModel does not have yet (capacity basal BC, 2nd-order thermal advection, ...) take
    # the ported choices, logged by with_ported_options: the Julia backend simulates that variant.
    return P.with_ported_options(P.YelmoParameters("Greenland", (getfield(p0, g) for g in P.GROUPS)...))
end

function setup_warm()
    p = ismip7_parameters_julia()
    # N_eff set externally (by the hydrology and its u_b hook): Yelmo.jl's yhyd.bkt_N_closure = -1
    p = _override_field(p, :yhyd, _override_field(p.yhyd, :bkt_N_closure, -1))
    mkpath(RUN_DIR)
    yelmo = YelmoModel(ISMIP7_RESTART, 0.0; p, boundaries = :bounded, rundir = RUN_DIR, strict = false)
    @info "warm start" restart = ISMIP7_RESTART ytill_method = p.ytill.method solver = p.ydyn.solver
    Yelmo.update_diagnostics!(yelmo)
    Yelmo.YelmoModelDyn.dyn_step!(yelmo, 0.0)
    return yelmo
end

function main_aguamenti(outdir, dt, time_end)
    mkpath(outdir)
    yelmo = get(ENV, "WARM", "1") == "1" ? setup_warm() : setup_yelmo(; external_neff = true)
    sim = HYDRO_MODEL == "aguamenti" ? build_hydrology_sim_Aguamenti(yelmo) : build_hydrology_sim_K24(yelmo)
    coupling = CoupledHydrology(sim)
    install_ub_hook!(coupling, yelmo)
    agu = sim.model isa AguamentiHydroModel
    @info "coupled run" HYDRO_MODEL hook = (yelmo.hooks.neff_from_ub !== nothing) dt time_end AGU_KAPPA AGU_COLD_T threads = Threads.nthreads()

    nx, ny = yelmo.g.Nx, yelmo.g.Ny
    nsteps = round(Int, time_end / dt)
    ds = NCDataset(joinpath(outdir, "yelmo.nc"), "c")
    defDim(ds, "x", nx); defDim(ds, "y", ny); defDim(ds, "time", nsteps + 1)
    vt = defVar(ds, "time", Float64, ("time",))
    v = Dict(n => defVar(ds, n, Float64, ("x", "y", "time")) for n in
             (agu ? ("N_eff", "f_grnd", "uxy_b", "taub", "T_prime_b", "H_w", "regime", "chi", "N_t", "dry", "q") :
                    ("N_eff", "f_grnd", "uxy_b", "taub", "T_prime_b", "H_w", "q")))
    vs = Dict(n => defVar(ds, n, Float64, ("time",)) for n in ("step_s", "hydro_s", "ssa_iters", "hook_calls", "hydro_iters", "hydro_conv"))
    hydro_s = Ref(0.0)
    function save(k, step_s, ssa)
        vt[k] = yelmo.time
        v["N_eff"][:, :, k]  = interior(yelmo.dyn.N_eff, :, :, 1)
        v["f_grnd"][:, :, k] = interior(yelmo.tpo.f_grnd, :, :, 1)
        v["uxy_b"][:, :, k]  = interior(yelmo.dyn.uxy_b, :, :, 1)
        v["taub"][:, :, k]   = interior(yelmo.dyn.taub, :, :, 1)
        v["T_prime_b"][:, :, k] = interior(yelmo.thrm.T_prime_b, :, :, 1)
        v["H_w"][:, :, k]    = interior(yelmo.thrm.H_w, :, :, 1)
        v["q"][:, :, k]      = Array(sim.model.q)
        if agu
            v["regime"][:, :, k] = sim.model.regime
            v["chi"][:, :, k]    = sim.model.chi
            v["N_t"][:, :, k]    = sim.model.N_t
            v["dry"][:, :, k]    = sim.model.dry
        end
        vs["step_s"][k] = step_s; vs["hydro_s"][k] = hydro_s[]; vs["ssa_iters"][k] = ssa
        vs["hook_calls"][k] = NEFF_CALLS[]
        vs["hydro_iters"][k] = agu ? sim.model.info.iters : -1
        vs["hydro_conv"][k]  = agu ? Float64(sim.model.info.converged) : -1
    end
    save(1, 0.0, 0)
    for k in 1:nsteps
        t0 = time(); NEFF_CALLS[] = 0
        # the hydrology update and Yelmo's step, timed separately
        hydro_s[] = @elapsed couple_step!(sim.model, sim, yelmo, dt)
        Yelmo.step!(yelmo, dt)
        ssa = hasproperty(yelmo.dyn, :scratch) ? yelmo.dyn.scratch.ssa_iter_now[] : -1
        save(k + 1, time() - t0, ssa)
        g = interior(yelmo.tpo.f_grnd, :, :, 1) .> 0.5
        if agu
            m = sim.model
            r = round.(Int, m.regime[g]); n = length(r)
            frac(c) = round(count(==(c), r) / n; digits = 3)
            @info "step $k t=$(round(yelmo.time; digits = 3)) yr  $(round(time() - t0; digits = 1)) s (hydro $(round(hydro_s[]; digits = 2)) s, $(m.info.iters) sweeps, $(m.info.converged ? "ok" : "NOT converged"))  SSA iters $ssa  hook calls $(NEFF_CALLS[])  " *
                  "N mean $(round(mean(interior(yelmo.dyn.N_eff, :, :, 1)[g]) / 1e6; digits = 3)) MPa  cold $(frac(0)) dry $(frac(1)) uplift $(frac(2)) press $(frac(3)) free $(frac(4)) strength $(frac(5))  chi>0.5 $(count(m.chi[g] .> 0.5))"
        else
            @info "step $k t=$(round(yelmo.time; digits = 3)) yr  $(round(time() - t0; digits = 1)) s (hydro $(round(hydro_s[]; digits = 2)) s)  SSA iters $ssa  hook calls $(NEFF_CALLS[])  N mean $(round(mean(interior(yelmo.dyn.N_eff, :, :, 1)[g]) / 1e6; digits = 3)) MPa"
        end
        all(isfinite, interior(yelmo.dyn.uxy_b, :, :, 1)) || error("non-finite basal velocity at step $k")
    end
    close(ds)
end

main_aguamenti(ARGS[1], parse(Float64, get(ARGS, 2, "0.1")), parse(Float64, get(ARGS, 3, "2.0")))
