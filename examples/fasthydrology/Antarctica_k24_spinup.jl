## ANT-32KM ISMIP7 spin-up with the K24 subglacial hydrology (Kazmierczak et al. 2024) from FastHydrology.jl
## driving YelmoMirror, reproducing the yelmox ISMIP7 spin-up j32_ours (Javi's parameters, 30 kyr,
## friction + basin-mean thermal-forcing optimisation) with N_eff from K24 instead of Yelmo's own closure.
##
## The loop is yelmox_esm's (yelmox main, yelmox_esm/yelmox_esm.f90), with the same compiled physics for
## everything but the hydrology and the optimisers: Yelmo through libyelmo_c_api.so, isostasy + sea level,
## ESM climate / direct SMB and marine melt through the yelmox C APIs (yelmox_capi.jl). The cb_ref and
## tf_corr optimisers are the Julia ports (Yelmo.optimize_cb_ref!, optimize_tf_corr_basin! below). One
## outer step of dtt (10 yr), as yelmox:
##
##   step_optimize -> [K24 update: push Yelmo fields, route, push N_eff, C_frz, Q_diss] -> isostasy ->
##   couple smb / marine -> Yelmo step -> climate + marine forcing
##
## Within the Yelmo step, Yelmo's DIVA iteration calls back into K24 (yelmo_set_neff_callback!) for N at
## each iteration's basal speed, routing held from the step's K24 update (install_ub_hook!).
##
## Start: the t = 0 restart bundle of the reference run, i.e. its cold start (present-day geometry,
## robin-cold temperature, cb_ref = cf_init after its t = 0 optimiser call, shelves killed, isostasy and sea
## level), entered through the classic restart branch. The validated route: the same driver without
## hydrology reproduced a 20 kyr yelmox spin-up to mean |dH| < 1 m (examples/ismip7_validation).
##
## Namelist: the reference run's own, verbatim, except
##   - yelmo.restart: the reference t = 0 bundle
##   - keys renamed since the reference run's yelmo, same meaning: ydyn.slide_T = F (was scale_T = 0),
##     ycalv.tau_ice_flt = tau_ice_grnd = 125e3 (was tau_ice), ytherm.strain_heating = "full" (was
##     use_strain_sia = F)
##   - ytherm: method = "enth", basal_bc_method = "capacity", qb_method = 1 (frictional heat on the C-grid
##     faces), cap_source = "hyd" (K24's freeze-on capacity) or "till" (control)
##   - yhyd.bkt_N_closure = -1 (N_eff from K24) for SPIN_HYDRO = k24
## and input/ holds the yelmo build's own defaults / variable tables / constants (Earth sec_year set to the
## reference run's 365 d), the rest of the reference run's input/ linked in.
##
## Environment:
##   SPIN_RUNDIR      run directory (required)
##   SPIN_HYDRO       k24 (default) | bucket (control: Yelmo's own N closure, same everything else)
##   SPIN_REF_DIR     reference yelmox run dir [~/yelmox/output/ismip7_ant/j32_ours]
##   SPIN_T_END       [30000] yr (model time);  SPIN_RESTART_DT [1000] yr;  SPIN_SNAP_DT [1000] yr
##   SPIN_OPT_T_END   [nml cf/tf_time_end] yr: keep the cb_ref / tf_corr optimisers on until this time
##   SPIN_START       restart directory of an earlier spin-up (write_restart: yelmo_restart.nc + tf_corr.nc)
##                    to continue from, at its own time; "" (default): the reference t = 0 bundle. The
##                    optimisers follow the reference schedule (&opt *_time_end), so a continuation past
##                    30 kyr is a free run with cb_ref and tf_corr held. Needs SPIN_ISOS=0.
##   SPIN_ISOS        1 (default): isostasy + sea level as yelmox; 0: bed and sea level held at the start
##                    state (the restarts written here hold no isostasy state)
##   YELMO_INPUT_DIR  input/ of the yelmo build libyelmo_c_api.so came from
##   YELMOX_CAPI_DIR  see yelmox_capi.jl
## plus the K24 knobs of Antarctica_yelmo-fasthydrology.jl (K24_KAPPA, K24_FRICTION, K24_SIGMAT, ...);
## FASTHYDRO_BACKEND must be "mirror".

ENV["KRYO_NO_MAIN"] = "1"
include(joinpath(@__DIR__, "Antarctica_yelmo-fasthydrology.jl"))   # K24 build / push / pull, DIVA callback, knobs
include(joinpath(@__DIR__, "yelmox_capi.jl"))
using Yelmo.YelmoMirrorPar: parse_nml_file
using Printf

BACKEND == "mirror" || error("Antarctica_k24_spinup.jl runs on YelmoMirror: set FASTHYDRO_BACKEND=mirror")

const SPIN_HYDRO  = lowercase(get(ENV, "SPIN_HYDRO", "k24"))
SPIN_HYDRO in ("k24", "bucket") || error("SPIN_HYDRO must be 'k24' or 'bucket', got '$SPIN_HYDRO'")
const REF_DIR     = get(ENV, "SPIN_REF_DIR", joinpath(homedir(), "yelmox", "output", "ismip7_ant", "j32_ours"))
const REF_NML     = joinpath(REF_DIR, "yelmox_esm_Antarctica_ismip7.nml")
const REF_BUNDLE  = joinpath(REF_DIR, "restart-0.000-kyr")
const SPIN_RUNDIR = ENV["SPIN_RUNDIR"]
const T_END       = parse(Float64, get(ENV, "SPIN_T_END", "30000"))
const RESTART_DT  = parse(Float64, get(ENV, "SPIN_RESTART_DT", "1000"))
const SNAP_DT     = parse(Float64, get(ENV, "SPIN_SNAP_DT", "1000"))
const YELMO_INPUT = get(ENV, "YELMO_INPUT_DIR", joinpath(homedir(), "wt", "yelmo-k24spin", "input"))
const SEC_YEAR_REF = 31536000.0   # the reference run's Earth sec_year (365 d)
# Initial cb_ref: "bundle" (default) the reference start's own; "match" rescaled by N_ref/N_K24 so that
# c_bed = cb_ref*N at t = 0 equals the reference start's (K24 N is far below the reference closure's);
# "scale<f>" multiplied by f. Clamped to [ytill.cf_min, ytill.cf_ref]. SPIN_CF_REF overrides ytill.cf_ref.
const SPIN_CB_INIT   = get(ENV, "SPIN_CB_INIT", "bundle")
const SPIN_CF_REF    = get(ENV, "SPIN_CF_REF", "")
const SPIN_THERM     = get(ENV, "SPIN_THERM", "enth")   # enth: enthalpy + capacity basal BC (default); temp: the reference run's ytherm as is (diagnostics)
const SPIN_QB_METHOD = parse(Int, get(ENV, "SPIN_QB_METHOD", "1"))   # Yelmo frictional heat; 1 = faces (pairs with K24_FRICTION=faces)
const SPIN_START     = get(ENV, "SPIN_START", "")
const SPIN_ISOS      = get(ENV, "SPIN_ISOS", "1") != "0"
isempty(SPIN_START) || SPIN_ISOS == false || error("SPIN_START needs SPIN_ISOS=0 (no isostasy state in the restarts)")
isempty(SPIN_START) || SPIN_CB_INIT == "bundle" || error("SPIN_START continues the restart's own cb_ref: SPIN_CB_INIT must be bundle")

fld(f)      = Array{Float64}(interior(f)[:, :, 1])
setf!(f, A) = (interior(f)[:, :, 1] .= A; f)

# ── run directory and namelist ────────────────────────────────────────────────

"""Remove `key` from namelist group `&group` of `txt` (an error if it is not there)."""
function _nml_del(txt::AbstractString, group::AbstractString, key::AbstractString)
    m = match(Regex("^&$(group)[ \\t]*\\n(.*?)^/", "ms"), txt)
    m === nothing && error("namelist group &$group not found")
    lines = split(chomp(m.captures[1]), '\n')
    keep = filter(l -> (k = match(r"^\s*(\w+)\s*=", l); k === nothing || lowercase(k.captures[1]) != lowercase(key)), lines)
    length(keep) == length(lines) - 1 || error("&$group $key: expected exactly one entry")
    return txt[1:prevind(txt, m.offset)] * "&$(group)\n" * join(keep, "\n") * "\n/" * txt[m.offset + ncodeunits(m.match):end]
end

# Keys of the reference namelist the yelmo build no longer has (it stops on unknown keys). Renamed ones are
# set under their new name in spinup_nml; the rest belong to code that was removed.
const NML_OBSOLETE = [("yelmo", "cfl_diff_max"), ("yelmo", "pc_corr_vel"), ("ycalv", "tau_ice"),
    ("ydyn", "cb_sia"), ("ydyn", "scale_T"), ("ydyn", "ssa_beta_max"), ("ydyn", "T_frz"),
    ("ytherm", "cp_rock"), ("ytherm", "use_strain_sia"), ("ytopo", "dHdt_dyn_lim"), ("ytopo", "f_ice_method"),
    ("ytopo", "margin2nd"), ("ytopo", "margin_flt_subgrid"), ("ytopo", "surf_gl_method"),
    ("yhyd", "k24_long_coupling_water"),
    # Fortran-internal K24 (unused: N_eff comes from FastHydrology.jl or the bucket), but validated at init;
    # the reference value 0 conflicts with the build's default routing (Warner needs the taped solver, 3)
    ("yhyd", "k24_flux_solver")]

"""The reference run's namelist with the edits listed in the header."""
function spinup_nml()
    txt = read(REF_NML, String)
    for (g, k) in NML_OBSOLETE
        txt = _nml_del(txt, g, k)
    end
    txt = _nml_set(txt, "yelmo", "restart", joinpath(isempty(SPIN_START) ? REF_BUNDLE : SPIN_START, "yelmo_restart.nc"))
    txt = _nml_set(txt, "ydyn", "slide_T", false)
    txt = _nml_set(txt, "ycalv", "tau_ice_flt", 125e3)
    txt = _nml_set(txt, "ycalv", "tau_ice_grnd", 125e3)
    txt = _nml_set(txt, "ytherm", "strain_heating", "full")
    if SPIN_THERM == "enth"
        txt = _nml_set(txt, "ytherm", "method", "enth")
        txt = _nml_set(txt, "ytherm", "basal_bc_method", "capacity")
        txt = _nml_set(txt, "ytherm", "cap_source", SPIN_HYDRO == "k24" ? "hyd" : "till")
    end
    txt = _nml_set(txt, "ytherm", "qb_method", SPIN_QB_METHOD)
    SPIN_HYDRO == "k24" && (txt = _nml_set(txt, "yhyd", "bkt_N_closure", -1))
    isempty(SPIN_CF_REF) || (txt = _nml_set(txt, "ytill", "cf_ref", parse(Float64, SPIN_CF_REF)))
    return txt
end

_link!(src, dst) = (ispath(dst) || islink(dst)) || symlink(src, dst)

function prepare_rundir(rundir)
    mkpath(rundir)
    for d in ("ice_data", "isostasy_data", "maps")
        _link!(realpath(joinpath(REF_DIR, d)), joinpath(rundir, d))
    end
    inp = joinpath(rundir, "input"); mkpath(inp)
    own = readdir(YELMO_INPUT)
    for f in own
        f == "yelmo_phys_const.nml" || _link!(joinpath(YELMO_INPUT, f), joinpath(inp, f))
    end
    refinp = realpath(joinpath(REF_DIR, "input"))
    for f in readdir(refinp)
        (f in own || startswith(f, "yelmo")) && continue   # yelmo's own files come from the build
        src = joinpath(refinp, f)
        ispath(src) || continue                            # dangling link in the reference input/
        _link!(realpath(src), joinpath(inp, f))
    end
    pc = _nml_set(read(joinpath(YELMO_INPUT, "yelmo_phys_const.nml"), String), "Earth", "sec_year", SEC_YEAR_REF)
    write(joinpath(inp, "yelmo_phys_const.nml"), pc)
    write(joinpath(rundir, "run.nml"), spinup_nml())
    return joinpath(rundir, "run.nml")
end

# ── namelist access ───────────────────────────────────────────────────────────

nmlstr(nml, g, k)  = String(strip(replace(nml[g][k], "'" => "", "\"" => "")))
nmlflt(nml, g, k)  = parse(Float64, replace(nmlstr(nml, g, k), r"[dD]" => "e"))
nmlbool(nml, g, k) = occursin(r"^(\.?true\.?|t)$"i, nmlstr(nml, g, k))
nmlints(nml, g, k) = parse.(Int, split(replace(nmlstr(nml, g, k), "," => " ")))

# ── optimisers (yelmox step_optimize) ─────────────────────────────────────────

struct OptPars
    opt_cf::Bool; cf_t0::Float64; cf_t1::Float64; tau_c::Float64; H0::Float64; fill_method::String
    cf_min::Float64; cf_max::Float64
    opt_tf::Bool; tf_t0::Float64; tf_t1::Float64; H_grnd_lim::Float64; tau_m::Float64; m_temp::Float64
    tf_min::Float64; tf_max::Float64; tf_basins::Vector{Int}
    rel_time2::Float64
end

# SPIN_OPT_T_END overrides the nml cf/tf_time_end, e.g. to keep optimising in a continuation past 30 kyr
_opt_t_end(t_nml) = haskey(ENV, "SPIN_OPT_T_END") ? parse(Float64, ENV["SPIN_OPT_T_END"]) : t_nml

# yelmox domain_opt_init: the cb_ref bounds are Yelmo's till_cf_min / till_cf_ref
load_opt(nml) = OptPars(
    nmlbool(nml, "opt", "opt_cf"), nmlflt(nml, "opt", "cf_time_init"), _opt_t_end(nmlflt(nml, "opt", "cf_time_end")),
    nmlflt(nml, "opt", "tau_c"), nmlflt(nml, "opt", "H0"), nmlstr(nml, "opt", "fill_method"),
    nmlflt(nml, "ytill", "cf_min"), nmlflt(nml, "ytill", "cf_ref"),
    nmlbool(nml, "opt", "opt_tf"), nmlflt(nml, "opt", "tf_time_init"), _opt_t_end(nmlflt(nml, "opt", "tf_time_end")),
    nmlflt(nml, "opt", "H_grnd_lim"), nmlflt(nml, "opt", "tau_m"), nmlflt(nml, "opt", "m_temp"),
    nmlflt(nml, "opt", "tf_min"), nmlflt(nml, "opt", "tf_max"), nmlints(nml, "opt", "tf_basins"),
    nmlflt(nml, "opt", "rel_time2"))

"""Basin-mean thermal-forcing correction: port of `optimize_tf_corr_basin` (yelmo libs/ice_optimization.f90),
the update yelmox uses with YELMOX_TF_BASIN=1. Per basin, over the cells near flotation
(H_grnd < H_grnd_lim) with observed or modelled ice, the mean thickness error and dH/dt give one rate
applied to the whole basin. `tf_basins[1] < 0`: every basin id present in `basins`."""
function optimize_tf_corr_basin!(tf_corr, H_ice, H_grnd, dHdt, H_obs, basins, H_grnd_lim, tau_m, m_temp,
                                 tf_min, tf_max, tf_basins, dt)
    f_damp = 2.0
    tol    = 1e-5
    blist  = tf_basins[1] < 0 ? sort(unique(vec(basins))) : Float64.(filter(>(0), tf_basins))
    for b in blist
        mask = (abs.(basins .- b) .< tol) .& (H_grnd .< H_grnd_lim) .& ((H_obs .> 0) .| (H_ice .> 0))
        n = count(mask)
        n > 0 || continue
        H_err = sum(H_ice[mask]) / n - sum(H_obs[mask]) / n
        dHdt_bar = sum(dHdt[mask]) / n
        rate = 1 / (tau_m * m_temp) * (H_err / tau_m + f_damp * dHdt_bar)
        tf_corr[basins .== b] .+= rate * dt
    end
    clamp!(tf_corr, tf_min, tf_max)
    return tf_corr
end

function step_optimize!(y, o::OptPars, elapsed, dt, dx)
    # topography relaxation: yelmox only applies it while elapsed <= rel_time2 (1 yr here), i.e. never
    # after the t = 0 step that is already in the start bundle; the nml's topo_rel = 0 holds from then on
    elapsed > o.rel_time2 || error("restart before rel_time2 is not supported")
    H_ice = fld(y.tpo.H_ice); dHdt = fld(y.tpo.dHidt); H_obs = fld(y.dta.pd_H_ice)
    if o.opt_cf && o.cf_t0 <= elapsed <= o.cf_t1
        cb = fld(y.dyn.cb_ref)
        optimize_cb_ref!(cb, H_ice, dHdt, fld(y.dyn.ux_s), fld(y.dyn.uy_s), H_obs, fld(y.dta.pd_uxy_s),
                         fld(y.dta.pd_H_grnd), o.cf_min, o.cf_max, dx, o.tau_c, o.H0, dt;
                         fill_method = o.fill_method, cb_tgt = fld(y.dyn.cb_tgt))
        setf!(y.dyn.cb_ref, cb)
    end
    if o.opt_tf && o.tf_t0 <= elapsed <= o.tf_t1
        tf = marshelf_get_var2D!(zeros(size(H_ice)), "tf_corr")
        optimize_tf_corr_basin!(tf, H_ice, fld(y.tpo.H_grnd), dHdt, H_obs, fld(y.bnd.basins), o.H_grnd_lim,
                                o.tau_m, o.m_temp, o.tf_min, o.tf_max, o.tf_basins, dt)
        marshelf_set_var2D!(tf, "tf_corr")
    end
    return nothing
end

"""Present-day bed as Yelmo loads it (&yelmo_data pd_topo_path, second of pd_topo_names)."""
function load_pd_bed(nml, domain, grid_name)
    path = replace(nmlstr(nml, "yelmo_data", "pd_topo_path"), "{domain}" => domain, "{grid_name}" => grid_name)
    name = strip.(split(replace(nml["yelmo_data"]["pd_topo_names"], "'" => "", "\"" => "")))[2]
    return NCDataset(ds -> Array{Float64}(coalesce.(ds[name][:, :], NaN)), path)
end

# ── forcing and the Yelmo step (yelmox couplers) ──────────────────────────────

function climate_and_marine!(y, t, dx)
    H = fld(y.tpo.H_ice); zb = fld(y.bnd.z_bed); fg = fld(y.tpo.f_grnd); zsl = fld(y.bnd.z_sl); bas = fld(y.bnd.basins)
    esm_step_climate!(t, fld(y.tpo.z_srf), H, zb, fg, zsl, bas, fld(y.dta.pd_z_srf))
    esm_step_marine!(H, zb, fg, zsl, fld(y.bnd.regions), bas, dx)
    return nothing
end

function couple_and_step!(y, dt, conv_we_ie, time_rel)
    nx, ny = size(fld(y.tpo.H_ice))
    t = y.time + dt
    if SPIN_ISOS
        dzcorr = zeros(nx, ny); Yelmo.yelmo_get_var2D!(dzcorr, Vector{UInt8}("bnd_dzbdt_corr\0"), y.calias)
        isos_update!(fld(y.tpo.H_ice), dzcorr, t, time_rel)
        setf!(y.bnd.z_bed, isos_get_var2D!(zeros(nx, ny), "z_bed"))
        setf!(y.bnd.z_sl,  isos_get_var2D!(zeros(nx, ny), "z_ss"))
    end
    setf!(y.bnd.smb_ref, esm_get_var2D!(zeros(nx, ny), "smb") .* conv_we_ie .* 1e-3)   # mm w.e./yr -> m i.e./yr
    setf!(y.bnd.T_srf,   esm_get_var2D!(zeros(nx, ny), "tsrf"))
    setf!(y.bnd.bmb_shlf, marshelf_get_var2D!(zeros(nx, ny), "bmb_shlf"))
    setf!(y.bnd.T_shlf,   marshelf_get_var2D!(zeros(nx, ny), "T_shlf"))
    Yelmo.step!(y, dt)
    return nothing
end

# ── diagnostics and output ────────────────────────────────────────────────────

const TS_COLS = ["time", "wall_s", "V_ice_km3", "A_ice_km2", "rmse_H_m", "rmse_logu", "n_grnd", "n_cb_cap", "n_cb_floor",
                 "cb_med", "N_med_MPa", "N_p10_MPa", "N_p90_MPa", "n_N_le0", "frac_temperate", "mean_abs_dHdt", "max_uxy_s",
                 "tf_corr_mean", "neff_calls",
                 "s_opt", "s_k24", "s_yelmo", "s_hook", "s_forcing", "s_diag"]   # wall time of this outer step's parts [s]; s_hook is inside s_yelmo

function diagnostics(y, o, dx, wall)
    H = fld(y.tpo.H_ice); Hobs = fld(y.dta.pd_H_ice); g = (fld(y.tpo.f_grnd) .> 0.5) .& (H .> 0)
    ice = (H .> 0) .| (Hobs .> 0)
    u = fld(y.dyn.uxy_s); uo = fld(y.dta.pd_uxy_s); uv = (H .> 0) .& (uo .> 0) .& (u .> 0)
    cb = fld(y.dyn.cb_ref)[g]; N = fld(y.dyn.N_eff)[g] ./ 1e6
    q(v, p) = isempty(v) ? NaN : quantile(v, p)
    tf = marshelf_get_var2D!(zeros(size(H)), "tf_corr")
    return [y.time, wall, sum(H) * dx^2 * 1e-9, count(H .> 0) * dx^2 * 1e-6, sqrt(mean((H .- Hobs)[ice] .^ 2)),
            sqrt(mean((log10.(u[uv]) .- log10.(uo[uv])) .^ 2)), count(g), count(cb .>= 0.999 * o.cf_max),
            count(cb .<= 1.001 * o.cf_min), q(cb, 0.5), q(N, 0.5), q(N, 0.1), q(N, 0.9), count(N .<= 0),
            mean(fld(y.thrm.T_prime_b)[g] .>= -0.1), mean(abs.(fld(y.tpo.dHidt)[H .> 0])), maximum(u),
            mean(tf[H .> 0]), NEFF_CALLS[]]
end

const SNAP_VARS = [("H_ice", y -> y.tpo.H_ice), ("z_srf", y -> y.tpo.z_srf), ("z_bed", y -> y.bnd.z_bed),
                   ("f_grnd", y -> y.tpo.f_grnd), ("dHidt", y -> y.tpo.dHidt), ("uxy_s", y -> y.dyn.uxy_s),
                   ("uxy_b", y -> y.dyn.uxy_b), ("taub", y -> y.dyn.taub), ("N_eff", y -> y.dyn.N_eff),
                   ("cb_ref", y -> y.dyn.cb_ref), ("c_bed", y -> y.dyn.c_bed), ("beta", y -> y.dyn.beta),
                   ("T_prime_b", y -> y.thrm.T_prime_b), ("bmb_grnd", y -> y.thrm.bmb_grnd), ("Q_b", y -> y.thrm.Q_b),
                   ("Q_ice_b", y -> y.thrm.Q_ice_b), ("smb_ref", y -> y.bnd.smb_ref), ("T_srf", y -> y.bnd.T_srf),
                   ("bmb_shlf", y -> y.bnd.bmb_shlf), ("z_sl", y -> y.bnd.z_sl)]

function snapshot!(fn, y, coupling)
    nx, ny = size(fld(y.tpo.H_ice))
    new = !isfile(fn)
    NCDataset(fn, new ? "c" : "a") do ds
        if new
            defDim(ds, "xc", nx); defDim(ds, "yc", ny); defDim(ds, "time", Inf)
            defVar(ds, "time", Float64, ("time",))
            for (nm, _) in SNAP_VARS; defVar(ds, nm, Float32, ("xc", "yc", "time")); end
            for nm in ("tf_corr", "pd_H_ice", "pd_uxy_s")
                defVar(ds, nm, Float32, ("xc", "yc", "time"))
            end
            if coupling isa CoupledHydrology
                for nm in ("k24_W", "k24_q", "k24_mdot", "k24_Q_b", "k24_Q_diss", "k24_kappa", "k24_frozen")
                    defVar(ds, nm, Float32, ("xc", "yc", "time"))
                end
            end
        end
        k = length(ds["time"]) + 1
        ds["time"][k] = y.time
        for (nm, get) in SNAP_VARS; ds[nm][:, :, k] = fld(get(y)); end
        ds["tf_corr"][:, :, k]  = marshelf_get_var2D!(zeros(nx, ny), "tf_corr")
        ds["pd_H_ice"][:, :, k] = fld(y.dta.pd_H_ice)
        ds["pd_uxy_s"][:, :, k] = fld(y.dta.pd_uxy_s)
        if coupling isa CoupledHydrology
            m, st = coupling.sim.model, coupling.sim.state
            ds["k24_W"][:, :, k]      = interior(st.W, :, :, 1)
            ds["k24_q"][:, :, k]      = interior(m.q, :, :, 1)
            ds["k24_mdot"][:, :, k]   = interior(m.mdot_total, :, :, 1)
            ds["k24_Q_b"][:, :, k]    = interior(m.Q_b, :, :, 1)
            ds["k24_Q_diss"][:, :, k] = interior(m.Q_diss, :, :, 1)
            ds["k24_kappa"][:, :, k]  = interior(m.kappa, :, :, 1)
            ds["k24_frozen"][:, :, k] = K24_FROZEN[] === nothing ? zeros(nx, ny) : Float64.(K24_FROZEN[])
        end
    end
end

"""Yelmo restart + the marine-shelf tf_corr (the optimiser state yelmox keeps in marine_shelf.nc)."""
function write_restart(y, rundir)
    dir = joinpath(rundir, @sprintf("restart-%.3f-kyr", y.time / 1000)); mkpath(dir)
    yelmo_write_restart!(y, joinpath(dir, "yelmo_restart.nc"))
    tf = marshelf_get_var2D!(zeros(size(fld(y.tpo.H_ice))), "tf_corr")
    NCDataset(joinpath(dir, "tf_corr.nc"), "c") do ds
        defDim(ds, "xc", size(tf, 1)); defDim(ds, "yc", size(tf, 2))
        defVar(ds, "tf_corr", tf, ("xc", "yc"))
    end
    return dir
end

# ── main ──────────────────────────────────────────────────────────────────────

"""Run directory, YelmoMirror, couplers, K24 and the initial cb_ref: everything up to the time loop."""
function setup()
    nmlfile = prepare_rundir(SPIN_RUNDIR)
    cd(SPIN_RUNDIR)
    nml = parse_nml_file(nmlfile)
    o = load_opt(nml)
    @info "K24 spin-up" SPIN_HYDRO REF_DIR SPIN_RUNDIR T_END K24_RELAX_QT K24_TAU_N K24_W_SAT K24_COLD_ABSORB K24_FREEZE_DEMAND K24_FROZEN_BED K24_FROZEN_HYST SPIN_CB_INIT SPIN_CF_REF K24_KAPPA K24_KAPPA_BED K24_FRICTION K24_SIGMAT K24_KAMB86 K24_UB_HOOK K24_SLIDING K24_A_BASAL G_SOURCE opt = o

    p = YelmoMirrorParameters("k24spin")
    t0 = isempty(SPIN_START) ? 0.0 :
         NCDataset(ds -> Float64(ds["time"][end]), joinpath(SPIN_START, "yelmo_restart.nc"))
    y = YelmoMirror(p, t0; alias = "ylmo1", rundir = SPIN_RUNDIR, overwrite = true, nml_file = nmlfile)   # Fortran clock starts at t0
    init_state!(y, t0; thrm_method = "robin-cold")   # restart branch: the bundle's (restart's) state
    isapprox(y.time, t0; atol = 1e-6) || error("Yelmo time $(y.time) != start time $t0")

    H0 = fld(y.tpo.H_ice); nx, ny = size(H0)
    dx = abs(Float64(y.g.Δxᶜᵃᵃ))
    domain, grid_name = nmlstr(nml, "yelmo", "domain"), nmlstr(nml, "yelmo", "grid_name")
    dtt = nmlflt(nml, "spinup", "dtt")
    time_rel = nmlflt(nml, "spinup", "tstep_const") - 2000.0   # yelmox timeline_init(time_ref = 2000), "const"
    nmlbool(nml, "coupling", "lim_pd_ice") && error("lim_pd_ice = True is not wrapped")
    @info "grid" domain grid_name nx ny dx dtt time_rel

    # yelmox esm_cold_start: shelves beyond the present-day extent stay dead (already so in the bundle)
    if nmlbool(nml, "spinup", "kill_shelves")
        mk = fld(y.bnd.mask_ice); mk[fld(y.dta.pd_mask_bed) .== 0.0] .= 0.0; setf!(y.bnd.mask_ice, mk)
    end

    pcall = parse_nml_file(joinpath(SPIN_RUNDIR, "input", "yelmo_phys_const.nml"))
    pc = Dict(lowercase(k) => v for (k, v) in pcall[first(k for k in keys(pcall) if lowercase(k) == lowercase(nmlstr(nml, "yelmo", "phys_const")))])
    pcf(k) = parse(Float64, replace(strip(replace(pc[k], "'" => "")), r"[dD]" => "e"))
    conv_we_ie = pcf("rho_w") / pcf("rho_ice")

    # classic restart branch: marine shelf (tf_corr from the bundle), isostasy + sea level, ESM forcing
    marshelf_init!(nmlfile, "marine_shelf", nx, ny, domain, grid_name, fld(y.bnd.regions), fld(y.bnd.basins),
                   axis_centered_m(nx, dx), axis_centered_m(ny, dx), dx)
    tf0 = isempty(SPIN_START) ?
        NCDataset(ds -> Array{Float64}(coalesce.(ds["tf_corr"][:, :, 1], 0.0)), joinpath(REF_BUNDLE, "marine_shelf.nc")) :
        NCDataset(ds -> Array{Float64}(ds["tf_corr"][:, :]), joinpath(SPIN_START, "tf_corr.nc"))
    marshelf_set_var2D!(tf0, "tf_corr")
    esm_init!(nmlfile, SPIN_RUNDIR, domain, grid_name, nx, ny)
    if SPIN_ISOS
        isos_init!(nmlfile, "isos", nx, ny, dx, dx, time_rel)
        isos_init_state_restart!(REF_BUNDLE, fld(y.bnd.z_bed), H0, 0.0, time_rel)
        setf!(y.bnd.z_bed, isos_get_var2D!(zeros(nx, ny), "z_bed"))
        setf!(y.bnd.z_sl,  isos_get_var2D!(zeros(nx, ny), "z_ss"))
    end
    climate_and_marine!(y, t0, dx)
    maximum(abs.(esm_get_var2D!(zeros(nx, ny), "Qd_ann"))) == 0.0 || @warn "esm Qd_ann is non-zero (subglacial discharge is not coupled here)"

    if SPIN_HYDRO == "k24" && K24_KAPPA_BED == "pd"
        # Yelmo's present-day reference bed (bnd%z_bed_ref, the observed bed after Yelmo's topography
        # initialisation, before any isostatic change); cross-checked against the observation file
        K24_PD_BED[] = fld(y.bnd.z_bed_ref)
        gr = (fld(y.tpo.f_grnd) .> 0.5) .& (fld(y.tpo.H_ice) .> 0)
        dobs = abs.(load_pd_bed(nml, domain, grid_name)[gr] .- K24_PD_BED[][gr])
        median(dobs) < 1e-3 || error("z_bed_ref does not match the observed bed (median |diff| $(median(dobs)) m): grid mix-up?")
        κ = _kappa(y)
        @info "K24 bed type from the present-day bed (z_bed_ref)" cells_adjusted_by_yelmo_init = count(dobs .> 1) frac_hard_grounded = mean(κ[gr] .== 0) frac_soft_grounded = mean(κ[gr] .== 1) mean_kappa_grounded = mean(κ[gr])
    end
    coupling = SPIN_HYDRO == "k24" ? CoupledHydrology(build_hydrology_sim_K24(y)) : NoCoupling()
    install_ub_hook!(coupling, y)

    if SPIN_CB_INIT != "bundle"
        cb = fld(y.dyn.cb_ref); gr = (fld(y.tpo.f_grnd) .> 0.5) .& (fld(y.tpo.H_ice) .> 0)
        if SPIN_CB_INIT == "match"
            SPIN_HYDRO == "k24" || error("SPIN_CB_INIT=match needs SPIN_HYDRO=k24")
            N_ref = fld(y.dyn.N_eff)   # the reference closure's N in the start bundle
            couple_step!(coupling.sim.model, coupling.sim, y, dtt)   # K24 N on the start state (pushed again in step 1)
            N_k24 = fld(y.dyn.N_eff)
            fac = ifelse.(gr .& (N_k24 .> 0) .& (N_ref .> 0), N_ref ./ max.(N_k24, 1.0), 1.0)
        elseif startswith(SPIN_CB_INIT, "scale")
            fac = fill(parse(Float64, SPIN_CB_INIT[6:end]), size(cb))
        else
            error("SPIN_CB_INIT must be bundle, match or scale<f>, got $SPIN_CB_INIT")
        end
        cb_new = ifelse.(gr, clamp.(cb .* fac, o.cf_min, o.cf_max), cb)
        setf!(y.dyn.cb_ref, cb_new)
        @info "initial cb_ref rescaled" SPIN_CB_INIT median_factor = median(fac[gr]) p90_factor = quantile(fac[gr], 0.9) cb_med = median(cb_new[gr]) n_at_cap = count(cb_new[gr] .>= o.cf_max) n_at_floor = count(cb_new[gr] .<= o.cf_min)
    end

    return (; y, o, coupling, dx, dtt, conv_we_ie, time_rel, t0)
end

function main()
    (; y, o, coupling, dx, dtt, conv_we_ie, time_rel, t0) = setup()
    K24_RELAX_N_ON[] = true   # K24_TAU_N relaxation from the first loop step (the start-up N above is K24's own)
    tsfile = joinpath(SPIN_RUNDIR, "spinup_ts.tsv")
    open(tsfile, "w") do io; println(io, join(TS_COLS, '\t')); end
    snapfile = joinpath(SPIN_RUNDIR, "spinup_2D.nc")
    wall0 = time()
    row = vcat(diagnostics(y, o, dx, 0.0), zeros(6))
    open(tsfile, "a") do io; println(io, join(row, '\t')); end
    snapshot!(snapfile, y, coupling)

    nsteps = round(Int, (T_END - t0) / dtt)
    for n in 1:nsteps
        t = t0 + n * dtt
        hook0 = NEFF_TIME[]
        s_opt = @elapsed step_optimize!(y, o, t, dtt, dx)
        s_k24 = @elapsed (coupling isa CoupledHydrology && couple_step!(coupling.sim.model, coupling.sim, y, dtt))
        s_yelmo = @elapsed couple_and_step!(y, dtt, conv_we_ie, time_rel)   # isostasy + couplers + Yelmo step
        s_forcing = @elapsed climate_and_marine!(y, t, dx)
        s_hook = NEFF_TIME[] - hook0

        isapprox(y.time, t; atol = 1e-6) || error("Yelmo time $(y.time) != loop time $t")
        all(isfinite, fld(y.tpo.H_ice)) && all(isfinite, fld(y.dyn.uxy_s)) || error("non-finite state at t = $t")
        s_diag = @elapsed (row = diagnostics(y, o, dx, time() - wall0))
        row = vcat(row, [s_opt, s_k24, s_yelmo, s_hook, s_forcing, s_diag])
        open(tsfile, "a") do io; println(io, join(row, '\t')); end
        if n % 10 == 0 || n <= 5
            @printf("t=%7.0f  V=%.4e km3  rmseH=%6.1f m  rmse_logu=%.3f  cb@cap=%d  cb@floor=%d  N_med=%.2f MPa  wall=%.0f s\n",
                    row[1], row[3], row[5], row[6], row[8], row[9], row[11], row[2])
            flush(stdout)
        end
        abs(rem(t, SNAP_DT)) < 1e-6 && snapshot!(snapfile, y, coupling)
        abs(rem(t, RESTART_DT)) < 1e-6 && write_restart(y, SPIN_RUNDIR)
    end
    abs(rem(T_END, RESTART_DT)) < 1e-6 || write_restart(y, SPIN_RUNDIR)
    @info "spin-up finished" time = y.time wall_h = (time() - wall0) / 3600
    return y
end

get(ENV, "SPIN_NO_MAIN", "0") == "1" || main()
