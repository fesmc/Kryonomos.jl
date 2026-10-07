## ccall wrappers around the yelmox C APIs (fesmc/yelmox branch yelmox-c-api, libs/yelmox_*_c_api.f90):
## the classic yelmox_esm driver's own isostasy + sea level (FastIsostasy), ESM climate / direct SMB and
## ocean forcing (esm_forcing) and marine melt (marine_shelf), so a Julia driver around YelmoMirror runs
## the same compiled physics as yelmox. Build the libraries with libs/build_{isos,marshelf,esm}_c_api.sh
## after `make yelmox_esm` in a yelmox checkout whose nested yelmo is the one libyelmo_c_api.so came from
## (both link static libyelmo), and point YELMOX_CAPI_DIR at its libyelmox/include.
##
## Each library holds one module-level Fortran instance (no alias support), and assumes the identity-grid
## case: grid_isos == grid_clim == grid_mshlf == the Yelmo grid.

const YELMOX_CAPI_DIR = get(ENV, "YELMOX_CAPI_DIR", joinpath(homedir(), "wt", "yelmox-c-api", "libyelmox", "include"))
const isoslib     = joinpath(YELMOX_CAPI_DIR, "libyelmox_isos_c_api.so")
const esmlib      = joinpath(YELMOX_CAPI_DIR, "libyelmox_esm_c_api.so")   # marine_shelf + ESM forcing, one shared mshlf1
const marshelflib = esmlib

# ── isostasy + barystatic sea level ───────────────────────────────────────────

function isos_init!(filename::String, group::String, nx::Int, ny::Int, dx::Float64, dy::Float64, time_rel_init::Float64)
    ccall((:isos_init, isoslib), Cvoid, (Ptr{UInt8}, Ptr{UInt8}, Cint, Cint, Float64, Float64, Float64),
          filename * "\0", group * "\0", Cint(nx), Cint(ny), dx, dy, time_rel_init)
end

"""Restore isostasy + sea level from a yelmox restart bundle (the classic restart branch)."""
function isos_init_state_restart!(fldr::String, z_bed::Matrix{Float64}, H_ice::Matrix{Float64}, time::Float64, time_rel::Float64)
    nx, ny = size(z_bed)
    ccall((:isos_init_state_restart, isoslib), Cvoid, (Ptr{UInt8}, Ptr{Float64}, Ptr{Float64}, Float64, Float64, Cint, Cint),
          fldr * "\0", z_bed, H_ice, time, time_rel, Cint(nx), Cint(ny))
end

"""One isostasy step (bsl_update(time_rel) inside), loaded by `H_ice`, with Yelmo's `dzbdt_corr`."""
function isos_update!(H_ice::Matrix{Float64}, dwdt_corr::Matrix{Float64}, time::Float64, time_rel::Float64)
    nx, ny = size(H_ice)
    ccall((:isos_update, isoslib), Cvoid, (Ptr{Float64}, Ptr{Float64}, Float64, Float64, Cint, Cint),
          H_ice, dwdt_corr, time, time_rel, Cint(nx), Cint(ny))
end

function isos_get_var2D!(buffer::Matrix{Float64}, name::String)
    nx, ny = size(buffer)
    ccall((:isos_get_var2D, isoslib), Cvoid, (Ptr{Float64}, Cint, Cint, Ptr{UInt8}), buffer, Cint(nx), Cint(ny), name * "\0")
    return buffer
end

# ── marine_shelf ──────────────────────────────────────────────────────────────

function marshelf_init!(filename::String, group::String, nx::Int, ny::Int, domain::String, grid_name::String,
                        regions::Matrix{Float64}, basins::Matrix{Float64}, xc::Vector{Float64}, yc::Vector{Float64}, dx::Float64)
    ccall((:marshelf_init, marshelflib), Cvoid,
          (Ptr{UInt8}, Ptr{UInt8}, Cint, Cint, Ptr{UInt8}, Ptr{UInt8}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
          filename * "\0", group * "\0", Cint(nx), Cint(ny), domain * "\0", grid_name * "\0", regions, basins, xc, yc, dx)
end

function marshelf_get_var2D!(buffer::Matrix{Float64}, name::String)
    nx, ny = size(buffer)
    ccall((:marshelf_get_var2D, marshelflib), Cvoid, (Ptr{Float64}, Cint, Cint, Ptr{UInt8}), buffer, Cint(nx), Cint(ny), name * "\0")
    return buffer
end

function marshelf_set_var2D!(buffer::Matrix{Float64}, name::String)
    nx, ny = size(buffer)
    ccall((:marshelf_set_var2D, marshelflib), Cvoid, (Ptr{Float64}, Cint, Cint, Ptr{UInt8}), buffer, Cint(nx), Cint(ny), name * "\0")
    return buffer
end

"""Cell-centre axis of the ISMIP7 polar-stereographic grids [m] (ANT-32KM: nx = 191, -3040..3040 km)."""
axis_centered_m(n::Int, dx::Float64) = Float64[(i - 1 - (n - 1) / 2) * dx for i in 1:n]

# ── ESM climate / direct SMB and ocean forcing ────────────────────────────────

"""Read &ctrl/&esm/&<run_step> from `path_par` and initialise esm_forcing (direct-SMB path only)."""
function esm_init!(path_par::String, workdir::String, domain::String, grid_name::String, nx::Int, ny::Int)
    ccall((:esm_init, esmlib), Cvoid, (Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Cint, Cint),
          path_par * "\0", workdir * "\0", domain * "\0", grid_name * "\0", Cint(nx), Cint(ny))
end

"""yelmox step_climate_esm (identity grids): reference climatology, anomalies, direct SMB + T_srf."""
function esm_step_climate!(time::Float64, z_srf, H_ice, z_bed, f_grnd, z_sl, basins, pd_z_srf)
    nx, ny = size(H_ice)
    ccall((:esm_step_climate, esmlib), Cvoid,
          (Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Cint, Cint),
          time, z_srf, H_ice, z_bed, f_grnd, z_sl, basins, pd_z_srf, Cint(nx), Cint(ny))
end

"""yelmox step_marine_shelf_esm (identity grids): ocean at shelf depth + anomalies, then marshelf_update."""
function esm_step_marine!(H_ice, z_bed, f_grnd, z_sl, regions, basins, dx::Float64)
    nx, ny = size(H_ice)
    ccall((:esm_step_marine, esmlib), Cvoid,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64, Cint, Cint),
          H_ice, z_bed, f_grnd, z_sl, regions, basins, dx, Cint(nx), Cint(ny))
end

function esm_get_var2D!(buffer::Matrix{Float64}, name::String)
    nx, ny = size(buffer)
    ccall((:esm_get_var2D, esmlib), Cvoid, (Ptr{Float64}, Cint, Cint, Ptr{UInt8}), buffer, Cint(nx), Cint(ny), name * "\0")
    return buffer
end
