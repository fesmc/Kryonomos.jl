## Thin ccall wrapper around libyelmox_isos_c_api.so (fesmc/yelmox branch
## yelmox-c-api, libs/yelmox_isos_c_api.f90) -- the classic yelmox_esm
## driver's own isostasy physics (FastIsostasy's fastisostasy.f90 +
## barysealevel.f90, method=2/interactive_sealevel="quad", ANT-32KM's real
## GIA rheology file), not reachable through yelmo's own C API since it's
## yelmox-level code, not yelmo core. Same ccall conventions as
## marshelf_wrapper.jl / Yelmo.jl's own YelmoMirrorCoreFields.jl.
##
## Single global Fortran instance (isos1/bsl1, module-level in the .f90) --
## no alias support. Assumes grid_isos == grid_yelmo (true for the ANT-32KM
## ISMIP7 setup this was built for), so isos%out fields are used directly
## at the yelmo grid, no remap.

const isoslib = joinpath(homedir(), "yelmox-worktrees", "yelmox-c-api",
                          "libyelmox", "include", "libyelmox_isos_c_api.so")

function isos_init!(filename::String, group::String, nx::Int, ny::Int,
                     dx::Float64, dy::Float64, time_init::Float64)
    ccall((:isos_init, isoslib), Cvoid,
        (Ptr{UInt8}, Ptr{UInt8}, Cint, Cint, Float64, Float64, Float64),
        filename * "\0", group * "\0", Cint(nx), Cint(ny), dx, dy, time_init)
    return nothing
end

function isos_init_ref!(z_bed_ref::Matrix{Float64}, H_ice_ref::Matrix{Float64})
    nx, ny = size(z_bed_ref)
    ccall((:isos_init_ref, isoslib), Cvoid,
        (Ptr{Float64}, Ptr{Float64}, Cint, Cint),
        z_bed_ref, H_ice_ref, Cint(nx), Cint(ny))
    return nothing
end

function isos_init_state!(z_bed::Matrix{Float64}, H_ice::Matrix{Float64}, time::Float64)
    nx, ny = size(z_bed)
    ccall((:isos_init_state, isoslib), Cvoid,
        (Ptr{Float64}, Ptr{Float64}, Float64, Cint, Cint),
        z_bed, H_ice, time, Cint(nx), Cint(ny))
    return nothing
end

function isos_update!(H_ice::Matrix{Float64}, time::Float64)
    nx, ny = size(H_ice)
    ccall((:isos_update, isoslib), Cvoid,
        (Ptr{Float64}, Float64, Cint, Cint),
        H_ice, time, Cint(nx), Cint(ny))
    return nothing
end

function isos_get_var2D!(buffer::Matrix{Float64}, name::String)
    nx, ny = size(buffer)
    ccall((:isos_get_var2D, isoslib), Cvoid,
        (Ptr{Float64}, Cint, Cint, Ptr{UInt8}),
        buffer, Cint(nx), Cint(ny), name * "\0")
    return buffer
end
