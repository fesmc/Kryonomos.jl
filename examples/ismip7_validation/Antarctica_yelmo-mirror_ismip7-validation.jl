## YelmoMirror vs. the classic yelmox_esm driver -- no hydrology coupling yet.
#
# Goal: reproduce, from Julia, what the classic Fortran ISMIP7 spin-up
# (`yelmox_esm/run_ismip7_antarctica.sh spinup` in the `yelmox` repo, ANT-32KM)
# does -- same grid, same real RACMO forcing, same friction-coefficient
# optimizer (`Yelmo.optimize_cb_ref!`, ported from `ice_optimization.f90` in
# fesmc/yelmo PR #8 / Yelmo.jl PR #88) -- and compare the result against that
# run's restart bundle. This isolates whether the Mirror + optimizer port is
# trustworthy BEFORE any K24/Shakti hydrology gets layered on top (see
# `examples/fasthydrology/` for that, once this passes).
#
# STATUS (2026-09-24): runs end-to-end (real ANT-32KM grid, confirmed via
# YelmoMirrorParameters(filename, "Antarctica") -- the one-arg string
# YelmoMirror(filename, ...) constructor does NOT read the nml's domain/grid,
# it silently defaults to Greenland; this cost real debugging time). First
# comparison against the classic 200-yr reference showed a LARGE mismatch
# (Mirror H_ice mean 725.6 m vs. reference 1.0 m -- the classic run's own
# spin-up is deep into a collapse by t=200, see project memory), not
# numerical noise -- don't treat this script as validated until that's
# understood and a re-run comparison is closer.
#
# Two things still need verifying on a real run (marked VERIFY below):
#   1. Whether writing `interior(y.dyn.cb_ref) .= ...` on the Julia side and
#      then calling `step!` actually pushes it into Fortran before the step
#      (matching how `hyd_N` external-N coupling pushes via `yelmo_sync!`),
#      or whether an explicit `yelmo_sync!(y)` call is needed first.
#   2. The Southern-Hemisphere lapse-rate formula below is transcribed from
#      `esm_clim_update` in yelmox's `libs/esm_forcing.f90` -- worth diffing
#      against that source again once running, in case this file has since
#      changed relative to what was read while writing this script.
#
# Known simplifications vs. the classic run (flagged, not silently dropped):
#   - Isostasy: wired via libyelmox_isos_c_api.so (fesmc/yelmox branch
#     yelmox-c-api, libs/yelmox_isos_c_api.f90) -- the classic driver's own
#     FastIsostasy coupling (method=2, interactive_sealevel, real ANT-32KM
#     GIA rheology file), not a Julia reimplementation. Same rationale as
#     marine-melt below: FastIsostasy.jl (Kryonomos.jl's own
#     examples/fastisostasy/) is a DIFFERENT implementation of the same
#     algorithm and would not reproduce the classic run bit-for-bit, so this
#     wraps the real compiled physics instead.
#   - Marine melt / PICO (bmb_method="quad-nl" in the classic run): wired via
#     libyelmox_marshelf_c_api.so, same pattern -- see marshelf_wrapper.jl.
#   - Reads the SAME height-dimension-fixed RACMO file the classic run uses
#     (`~/yelmox/data_fixes/ANT-32KM_ERA5-3H_RACMO2.3p2_1979-2022_monthly_smb-fixed.nc`),
#     not a locally-copied one -- keep that file where it is, or update
#     CLIMATE_FILE below.

## Preamble ############################################################
cd(dirname(@__DIR__))
import Pkg; Pkg.activate(".")
#######################################################################

using Yelmo
using Yelmo: YelmoMirror

using NCDatasets
using Statistics
using Printf

# ── Paths ────────────────────────────────────────────────────────────────────
const RUN_DIR       = @__DIR__
const YELMO_NML_MIRROR = joinpath(RUN_DIR, "Antarctica_ismip7_mirror.nml")
const CLIMATE_FILE  = joinpath(homedir(), "yelmox", "data_fixes",
                                "ANT-32KM_ERA5-3H_RACMO2.3p2_1979-2022_monthly_smb-fixed.nc")
# The classic run's reference restart to compare against (see project memory
# for how it was produced -- `run_ismip7_antarctica.sh`-equivalent runme
# invocation, ANT-32KM, 200 yr, full colleague override set):
const REFERENCE_RESTART = joinpath(homedir(), "yelmox", "output", "ismip7_ant",
                                    "spinup_ANT-32KM", "restart-0.200-kyr", "yelmo_restart.nc")

# &marine_shelf params (bmb_method="quad-nl", gamma_quad_nl, etc.) -- reuse the classic
# run's own already-resolved nml rather than duplicating those values here; marshelf_init
# just needs a file with a &marine_shelf group in it, doesn't have to be the Mirror nml.
const MARSHELF_NML = joinpath(homedir(), "yelmox", "output", "ismip7_ant",
                               "spinup_ANT-32KM", "yelmox_esm_Antarctica_ismip7.nml")

include("marshelf_wrapper.jl")
include("isos_wrapper.jl")

# &isos / &barysealevel params (method=2, interactive_sealevel, ANT-32KM GIA
# rheology file) -- reuse the classic run's own already-resolved nml, same
# rationale as MARSHELF_NML above.
const ISOS_NML = MARSHELF_NML

# Reference climatology window, matching the classic nml's &spinup.time_ref = 1985.0, 2014.0.
# The fixed RACMO file's `time` axis is monthly, starting ~Jan 1979 (528 months total,
# 1979-2022 inclusive) -- index 72 = Jan 1985, index 431 = Dec 2014 (verified via ncks
# during setup: both land on physically sane ~262-267 K spatial-mean t2m).
const REF_TIME_START_IDX = 72    # Jan 1985 (0-based -> +1 for Julia's 1-based NCDatasets)
const REF_TIME_END_IDX   = 431   # Dec 2014

# ── Time stepping ────────────────────────────────────────────────────────────
const DT_OUTER_YR = 1.0     # matches the classic run's yelmo.dt_min=1.0
const T_END_YR    = 200.0   # matches the short validated classic run (spinup_ANT-32KM)

# ── Optimizer parameters (from the colleague's opt_ant.sh / the validated
#    runme override set -- see project_yelmo_hydro_optimizer.md) ─────────────
const OPT_TAU_C      = 100.0
const OPT_H0         = -1.0    # <=0 -> scaleH=true in optimize_cb_ref!
const OPT_CF_MIN     = 1e-4
const OPT_CF_MAX     = 1e0     # no explicit upper bound in the classic config; 1.0 is
                                 # effectively unconstrained for a friction coefficient here
const OPT_FILL_METHOD = "cf_min"

# tf_corr params, from the same validated runme override set (see project memory).
const OPT_TAU_M      = 10.0
const OPT_M_TEMP     = 10.0
const OPT_TF_MIN     = -1.0
const OPT_TF_MAX     = 1.0
const OPT_H_GRND_LIM = 500.0

# ── Southern-Hemisphere lapse rate (esm_clim_update, yelmox libs/esm_forcing.f90) ──
# lapse = lapse[1] + (lapse[2]-lapse[1]) * cos(2*pi*(month*30.4375 - 30.4375)/365.25)
# esm_ant_ismip7.nml: lapse = 8.0e-3 6.5e-3  [K/m, summer/winter halves as in the Fortran fmt]
const LAPSE = (8.0e-3, 6.5e-3)

function lapse_rate(month::Int)
    l1, l2 = LAPSE
    return l1 + (l2 - l1) * cos(2π * (month*30.4375 - 30.4375) / 365.25)
end

# ── Forcing: reference-period monthly climatology, elevation-lapse-corrected ──
"""
Reads t2m/smb/z_srf from CLIMATE_FILE over [REF_TIME_START_IDX, REF_TIME_END_IDX],
collapses to a 12-month climatology (plain mean over the 30-year window -- matches
`range_mean` in the Fortran `varslice_update` call for a monthly-periodic slice),
then lapse-corrects t2m to the model's actual surface elevation exactly as
`esm_clim_update` does. Returns (T_srf_2d, smb_2d) as annual means -- YelmoMirror's
bnd.T_srf/bnd.smb_ref are single 2D fields here (no monthly cycle kept in bnd),
consistent with how `apply_forcing_mirror!` in tutorial_yelmo_antarctica.jl works.
"""
function load_ismip7_reference_climatology(z_srf_model::AbstractMatrix)
    nx, ny = size(z_srf_model)
    t2m_clim = zeros(nx, ny, 12)
    smb_clim = zeros(nx, ny, 12)
    zs_ref   = zeros(nx, ny)

    NCDataset(CLIMATE_FILE) do ds
        t2m_all = ds["t2m"][:, :, (REF_TIME_START_IDX+1):(REF_TIME_END_IDX+1)]      # K
        smb_all = ds["smb"][:, :, (REF_TIME_START_IDX+1):(REF_TIME_END_IDX+1)]      # kg m^-2 month^-1
        zs_ref  .= Array{Float64}(coalesce.(ds["z_srf"][:, :], 0.0))                    # m, time-independent

        n_months = size(t2m_all, 3)
        for m in 1:12
            idx = m:12:n_months   # every m-th month across the 30-yr window
            t2m_clim[:, :, m] .= dropdims(mean(t2m_all[:, :, idx]; dims=3); dims=3)
            smb_clim[:, :, m] .= dropdims(mean(smb_all[:, :, idx]; dims=3); dims=3)
        end
    end

    # Elevation-lapse correct t2m to the model's actual surface (esm_clim_update):
    #   t2m(model) = t2m_ref + lapse(month) * (zs_ref - z_srf_model)
    T_srf_ann = zeros(nx, ny)
    smb_ann   = zeros(nx, ny)
    for m in 1:12
        lapse = lapse_rate(m)
        @. T_srf_ann += (t2m_clim[:, :, m] + lapse * (zs_ref - z_srf_model)) / 12
        @. smb_ann   += smb_clim[:, :, m] / 12   # no elevation correction for SMB (esm_clim_update: "No model elevation changes for SMB")
    end

    # smb units: kg m^-2 month^-1 -> m ice-equivalent/yr (matches &gcm_smb_ref's
    # `scaling = 12.0 0.0` in esm_ant_ismip7.nml: kg/m2/month -> mm/yr is *12, then
    # mm w.e./yr -> m i.e./yr needs /1000 * rho_w/rho_ice).
    RHO_W, RHO_ICE = 1000.0, 910.0
    smb_ann .= smb_ann .* 12.0 .* 1.0e-3 .* (RHO_W / RHO_ICE)

    return T_srf_ann, smb_ann
end

const OBS_FILE = joinpath(homedir(), "yelmox", "ice_data", "ISMIP7", "Antarctica", "ANT-32KM",
                          "obs", "ANT-32KM_ObsISMIP7-v1.1.nc")

"""
Real geothermal heat flux, matching the classic run's `ghf.obs_name =
geothermal_heat_flux1,geothermal_heat_flux2` (Staal et al. 2021 / a second product) --
averages the two products where both are valid, falls back to whichever one is valid
where only one is, matching the usual multi-product-blend convention for this field.
Both carry `_FillValue = -9e33`; anything that negative is treated as missing.
"""
function load_ghf()
    NCDataset(OBS_FILE) do ds
        g1 = Array{Float64}(replace(ds["geothermal_heat_flux1"][:, :], missing => NaN))
        g2 = Array{Float64}(replace(ds["geothermal_heat_flux2"][:, :], missing => NaN))
        valid1 = g1 .> -1.0e10
        valid2 = g2 .> -1.0e10
        ghf = zeros(size(g1))
        @. ghf = ifelse(valid1 & valid2, (g1 + g2) / 2,
                  ifelse(valid1, g1,
                  ifelse(valid2, g2, 55.0)))   # both missing -> same fallback constant as before
        return ghf
    end
end

function apply_forcing_mirror!(y)
    z_srf = interior(y.tpo.z_srf)[:, :, 1]
    T_srf, smb = load_ismip7_reference_climatology(z_srf)

    interior(y.bnd.T_srf)[:, :, 1]    .= T_srf
    interior(y.bnd.smb_ref)[:, :, 1]  .= smb
    interior(y.bnd.Q_geo)[:, :, 1]    .= load_ghf()

    fill!(interior(y.bnd.z_sl), 0.0)

    return y
end

# ── The optimizer driving loop ────────────────────────────────────────────────
"""
One outer step: advance the ice model, then nudge cb_ref toward matching the
observed present-day thickness/velocity (optimize_cb_ref!, ported from
ice_optimization.f90 -- see fesmc/Yelmo.jl PR #88 / project memory). Mirrors
how the classic driver only ever calls this from outside yelmo_update, once
per outer timestep -- not part of Yelmo.jl's own step! phase order.
"""
function step_with_optimizer!(y, dt, tf_corr::AbstractMatrix)
    step!(y, dt)   # VERIFY: confirm this both pushes any Julia-side cb_ref edits from
                    # the PREVIOUS iteration into Fortran before stepping, and pulls the
                    # post-step H_ice/dHidt/ux_bar/uy_bar back out afterward -- if not
                    # automatic, add an explicit yelmo_sync!(y) call here.

    H_ice   = interior(y.tpo.H_ice)[:, :, 1]
    dHdt    = interior(y.tpo.dHidt)[:, :, 1]
    ux_bar  = interior(y.dyn.ux_bar)[:, :, 1]
    uy_bar  = interior(y.dyn.uy_bar)[:, :, 1]
    H_obs      = interior(y.bnd.H_ice_ref)[:, :, 1]   # dta%pd%H_ice, copied to bnd%H_ice_ref at init
    uxy_obs    = interior(y.dta.pd_uxy_s)[:, :, 1]     # needs fesmc/yelmo PR #8's dta_pd_uxy_s getter
    H_grnd_obs = interior(y.dta.pd_H_grnd)[:, :, 1]    # needs fesmc/yelmo PR #8's dta_pd_H_grnd getter
    H_grnd     = interior(y.tpo.H_grnd)[:, :, 1]       # model-side (not obs) grounding state, for tf_corr
    cb_ref  = interior(y.dyn.cb_ref)[:, :, 1]

    optimize_cb_ref!(cb_ref, H_ice, dHdt, ux_bar, uy_bar, H_obs, uxy_obs, H_grnd_obs,
                      OPT_CF_MIN, OPT_CF_MAX, y.g.Δxᶜᵃᵃ, OPT_TAU_C, OPT_H0, dt;
                      fill_method=OPT_FILL_METHOD)

    interior(y.dyn.cb_ref)[:, :, 1] .= cb_ref   # write back; picked up by step! next iteration
                                                  # per the VERIFY note above

    # tf_corr now has a real consumer: libyelmox_marshelf_c_api.so wraps the classic
    # driver's own marine_shelf.f90 (fesmc/yelmox branch yelmox-c-api). Push tf_corr into
    # it, call the real compiled marshelf_update, pull the resulting bmb_shlf back out --
    # same physics the classic run uses, not a Julia reimplementation.
    optimize_tf_corr!(tf_corr, H_ice, H_grnd, dHdt, H_obs, OPT_H_GRND_LIM, y.g.Δxᶜᵃᵃ,
                       OPT_TAU_M, OPT_M_TEMP, OPT_TF_MIN, OPT_TF_MAX, dt)
    interior(y.bnd.T_shlf)[:, :, 1] .= tf_corr

    z_bed   = interior(y.bnd.z_bed)[:, :, 1]
    f_grnd  = interior(y.tpo.f_grnd)[:, :, 1]
    regions = interior(y.bnd.regions)[:, :, 1]
    basins  = interior(y.bnd.basins)[:, :, 1]
    z_sl    = interior(y.bnd.z_sl)[:, :, 1]
    dx      = abs(Float64(y.g.Δxᶜᵃᵃ))

    marshelf_set_var2D!(tf_corr, "tf_corr")
    marshelf_update!(H_ice, z_bed, f_grnd, regions, basins, z_sl, dx)

    bmb_shlf = similar(H_ice)
    marshelf_get_var2D!(bmb_shlf, "bmb_shlf")
    interior(y.bnd.bmb_shlf)[:, :, 1] .= bmb_shlf

    # Isostasy: step libyelmox_isos_c_api.so forward with the post-step H_ice,
    # pull the resulting bedrock/sea-surface fields back out, and write them
    # into bnd.z_bed/bnd.z_sl -- picked up by step! on the NEXT iteration,
    # exactly matching couple_isostasy_to_yelmo's landing point in the classic
    # driver (called before yelmo_update, i.e. before the next step).
    isos_update!(H_ice, y.time)
    z_bed_new = similar(H_ice)
    z_ss_new  = similar(H_ice)
    isos_get_var2D!(z_bed_new, "z_bed")
    isos_get_var2D!(z_ss_new,  "z_ss")
    interior(y.bnd.z_bed)[:, :, 1] .= z_bed_new
    interior(y.bnd.z_sl)[:, :, 1]  .= z_ss_new

    return y
end

function domain_diagnostics(y)
    H  = interior(y.tpo.H_ice)[:, :, 1]
    dx = abs(Float64(y.g.Δxᶜᵃᵃ))
    dy = abs(Float64(y.g.Δyᵃᶜᵃ))
    return (V_ice = sum(H) * dx * dy * 1e-9, H_max = maximum(H))
end

"""Compares this run's final H_ice against the classic yelmox_esm reference restart --
prints a summary; does not assert (this is exploratory parity-checking, not a test suite)."""
function compare_to_reference(y)
    isfile(REFERENCE_RESTART) || (@warn "Reference restart not found, skipping comparison" REFERENCE_RESTART; return)

    H_mirror = interior(y.tpo.H_ice)[:, :, 1]
    H_ref = NCDataset(REFERENCE_RESTART) do ds
        Array(ds["H_ice"][:, :])
    end

    if size(H_mirror) != size(H_ref)
        @warn "Grid size mismatch, cannot compare directly" size(H_mirror) size(H_ref)
        return
    end

    diff = H_mirror .- H_ref
    @printf("\n== Mirror vs. classic yelmox_esm reference (t=%.0f yr) ==\n", T_END_YR)
    @printf("  H_ice max   mirror=%.1f m   reference=%.1f m\n", maximum(H_mirror), maximum(H_ref))
    @printf("  H_ice mean  mirror=%.1f m   reference=%.1f m\n", mean(H_mirror), mean(H_ref))
    @printf("  |diff| mean=%.2f m   max=%.2f m   RMSE=%.2f m\n",
            mean(abs.(diff)), maximum(abs.(diff)), sqrt(mean(diff.^2)))
end

# ── Run ────────────────────────────────────────────────────────────────────────
function main()
    @info "Building Antarctica YelmoMirror (ISMIP7 validation, no hydrology)..."
    p = Yelmo.YelmoMirrorPar.YelmoMirrorParameters(YELMO_NML_MIRROR, "Antarctica")
    y = YelmoMirror(p, 0.0; alias="ismip7_ant_mirror", rundir=RUN_DIR, overwrite=true)

    apply_forcing_mirror!(y)
    init_state!(y, 0.0; thrm_method="robin-cold")

    nx, ny = size(interior(y.tpo.H_ice))[1:2]
    dx     = abs(Float64(y.g.Δxᶜᵃᵃ))
    regions = interior(y.bnd.regions)[:, :, 1]
    basins  = interior(y.bnd.basins)[:, :, 1]
    marshelf_init!(MARSHELF_NML, "marine_shelf", nx, ny, "Antarctica", "ANT-32KM",
                   regions, basins, axis_centered(nx, dx), axis_centered(ny, dx), dx)

    # Isostasy reference + initial state, mirroring yelmox_domain.f90's
    # domain_init_state: isos_init_ref (bed/thickness reference for the
    # elastic/viscous anomaly), THEN isos_init_state (current state at t=0).
    # z_bed_ref/H_ice_ref are real yelmo bnd fields (already exposed via the
    # C API as bnd_z_bed_ref/bnd_H_ice_ref, same generic getter H_ice_ref
    # already uses above for the cb_ref optimizer's H_obs).
    z_bed_ref = interior(y.bnd.z_bed_ref)[:, :, 1]
    H_ice_ref = interior(y.bnd.H_ice_ref)[:, :, 1]
    z_bed0    = interior(y.bnd.z_bed)[:, :, 1]
    H_ice0    = interior(y.tpo.H_ice)[:, :, 1]
    isos_init!(ISOS_NML, "isos", nx, ny, dx, dx, 0.0)
    isos_init_ref!(z_bed_ref, H_ice_ref)
    isos_init_state!(z_bed0, H_ice0, 0.0)

    d0 = domain_diagnostics(y)
    @info "t=0" V_ice_km3=d0.V_ice H_max_m=d0.H_max

    tf_corr = zeros(size(interior(y.tpo.H_ice))[1:2])   # caller-owned across iterations

    n_steps = round(Int, T_END_YR / DT_OUTER_YR)
    @printf("  %6s  %12s  %10s  %8s\n", "t[yr]", "V_ice[km³]", "H_max[m]", "step[s]")
    for k in 1:n_steps
        t_step = @elapsed step_with_optimizer!(y, DT_OUTER_YR, tf_corr)
        d = domain_diagnostics(y)
        @printf("  %6.0f  %12.4e  %10.1f  %8.2f\n", y.time, d.V_ice, d.H_max, t_step)
        flush(stdout)
    end

    compare_to_reference(y)

    return y
end

y = main()
