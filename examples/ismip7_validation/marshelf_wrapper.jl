## Thin ccall wrapper around libyelmox_marshelf_c_api.so (fesmc/yelmox branch
## yelmox-c-api, libs/yelmox_marshelf_c_api.f90) -- the classic yelmox_esm
## driver's own marine-melt physics (marine_shelf.f90, bmb_method="quad-nl"),
## not reachable through yelmo's own C API since it's yelmox-level code, not
## yelmo core. Same ccall conventions as Yelmo.jl's own
## YelmoMirrorCoreFields.jl (yelmo_get_var2D! etc) -- read that file first if
## this looks unfamiliar (Vector{UInt8}/String null-termination, Ptr{Float64}
## for arrays, Cint for dimensions).
##
## Single global Fortran instance (mshlf1, module-level in the .f90) -- no
## alias support, unlike YelmoMirror's ylmo1/ylmo2. Fine here: this script
## only ever drives one marine-shelf state per process.

const marshelflib = joinpath(homedir(), "yelmox-worktrees", "yelmox-c-api",
                              "libyelmox", "include", "libyelmox_marshelf_c_api.so")

function marshelf_init!(filename::String, group::String, nx::Int, ny::Int,
                         domain::String, grid_name::String,
                         regions::Matrix{Float64}, basins::Matrix{Float64},
                         xc::Vector{Float64}, yc::Vector{Float64}, dx::Float64)
    ccall((:marshelf_init, marshelflib), Cvoid,
        (Ptr{UInt8}, Ptr{UInt8}, Cint, Cint, Ptr{UInt8}, Ptr{UInt8},
         Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64),
        filename * "\0", group * "\0", Cint(nx), Cint(ny), domain * "\0", grid_name * "\0",
        regions, basins, xc, yc, dx)
    return nothing
end

function marshelf_update!(H_ice::Matrix{Float64}, z_bed::Matrix{Float64},
                           f_grnd::Matrix{Float64}, regions::Matrix{Float64},
                           basins::Matrix{Float64}, z_sl::Matrix{Float64}, dx::Float64)
    nx, ny = size(H_ice)
    ccall((:marshelf_update, marshelflib), Cvoid,
        (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
         Float64, Cint, Cint),
        H_ice, z_bed, f_grnd, regions, basins, z_sl, dx, Cint(nx), Cint(ny))
    return nothing
end

function marshelf_get_var2D!(buffer::Matrix{Float64}, name::String)
    nx, ny = size(buffer)
    ccall((:marshelf_get_var2D, marshelflib), Cvoid,
        (Ptr{Float64}, Cint, Cint, Ptr{UInt8}),
        buffer, Cint(nx), Cint(ny), name * "\0")
    return buffer
end

function marshelf_set_var2D!(buffer::Matrix{Float64}, name::String)
    nx, ny = size(buffer)
    ccall((:marshelf_set_var2D, marshelflib), Cvoid,
        (Ptr{Float64}, Cint, Cint, Ptr{UInt8}),
        buffer, Cint(nx), Cint(ny), name * "\0")
    return buffer
end

"""Symmetric cell-centre axis matching the ISMIP7 polar-stereographic grid
convention (verified against the classic run's log: nx=191, dx=32000 m ->
range -3040..3040 km). Yelmo's own grid object isn't threaded through to
this wrapper, so this is reconstructed from nx/dx directly rather than read
off `y.g`."""
axis_centered(n::Int, dx::Float64) = Float64[(i - 1 - (n - 1) / 2) * dx for i in 1:n]
