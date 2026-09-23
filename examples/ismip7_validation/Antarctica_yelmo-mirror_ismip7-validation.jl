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
# STATUS: first draft, not yet executed. Two things specifically need
# verifying on a real run (marked VERIFY below):
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
#   - Isostasy: OFF here (Antarctica_mirror.nml, which this is based on, has
#     no &isos section at all -- bed stays fixed). The classic run uses
#     LV-ELRA. Over 200 yr this is a second-order difference next to the
#     friction-optimization signal being tested, but matters for a longer run.
#   - tf_corr (thermal-forcing / ocean-melt optimization) is NOT driven here
#     yet -- only cb_ref. Add `Yelmo.optimize_tf_corr!` once cb_ref parity
#     is confirmed (see project_yelmo_hydro_optimizer.md memory).
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
using Yelmo.YelmoMirrorPar: YelmoParameters
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
        zs_ref  .= ds["z_srf"][:, :]                                                # m, time-independent

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

function apply_forcing_mirror!(y)
    z_srf = interior(y.tpo.z_srf)[:, :, 1]
    T_srf, smb = load_ismip7_reference_climatology(z_srf)

    interior(y.bnd.T_srf)[:, :, 1]    .= T_srf
    interior(y.bnd.smb_ref)[:, :, 1]  .= smb

    fill!(interior(y.bnd.Q_geo), 55.0)   # no per-cell GHF wired here yet; ISMIP7 obs file
                                          # has geothermal_heat_flux1/2 -- TODO: read those
                                          # instead, like the classic run's ghf.obs_name does.
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
function step_with_optimizer!(y, dt)
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
    cb_ref  = interior(y.dyn.cb_ref)[:, :, 1]

    optimize_cb_ref!(cb_ref, H_ice, dHdt, ux_bar, uy_bar, H_obs, uxy_obs, H_grnd_obs,
                      OPT_CF_MIN, OPT_CF_MAX, y.g.Δxᶜᵃᵃ, OPT_TAU_C, OPT_H0, dt;
                      fill_method=OPT_FILL_METHOD)

    interior(y.dyn.cb_ref)[:, :, 1] .= cb_ref   # write back; picked up by step! next iteration
                                                  # per the VERIFY note above

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
    p = YelmoParameters(YELMO_NML_MIRROR, "Antarctica")
    y = YelmoMirror(p, 0.0; alias="ismip7_ant_mirror", rundir=RUN_DIR, overwrite=true)

    apply_forcing_mirror!(y)
    init_state!(y, 0.0; thrm_method="robin-cold")

    d0 = domain_diagnostics(y)
    @info "t=0" V_ice_km3=d0.V_ice H_max_m=d0.H_max

    n_steps = round(Int, T_END_YR / DT_OUTER_YR)
    @printf("  %6s  %12s  %10s  %8s\n", "t[yr]", "V_ice[km³]", "H_max[m]", "step[s]")
    for k in 1:n_steps
        t_step = @elapsed step_with_optimizer!(y, DT_OUTER_YR)
        d = domain_diagnostics(y)
        @printf("  %6.0f  %12.4e  %10.1f  %8.2f\n", y.time, d.V_ice, d.H_max, t_step)
        flush(stdout)
    end

    compare_to_reference(y)

    return y
end

y = main()
