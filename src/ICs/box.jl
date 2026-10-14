export Box

"""
    Box{D, M, T} <: InitialCondition

Defines a discontinuous rectangular sub-region.

# Fields
- `u_bg::State{M, T}`: The background state applied outside the box.
- `u_box::State{M, T}`: The interior state applied within the box.
- `mins::Space{D, T}`: The lower spatial bounds of the rectangular region.
- `maxs::Space{D, T}`: The upper spatial bounds of the rectangular region.

Evaluates to `u_box` if `mins .<= pos .<= maxs` holds in every direction (the boundary
belongs to the box), and to `u_bg` otherwise. Mixed-precision inputs are promoted to one `T`.
"""
struct Box{D, M, T} <: InitialCondition
    u_bg::State{M, T}
    u_box::State{M, T}
    mins::Space{D, T}
    maxs::Space{D, T}
end

function Box(u_bg::State{M}, u_box::State{M}, mins::Space{D}, maxs::Space{D}) where {D, M}
    T = promote_type(eltype(u_bg), eltype(u_box), eltype(mins), eltype(maxs))
    return Box{D, M, T}(u_bg, u_box, mins, maxs)
end

(ic::Box)(pos::Space{D}) where {D} = all(ic.mins .<= pos .<= ic.maxs) ? ic.u_box : ic.u_bg

function build_ic(::Val{:box}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int
    
    u_bg = State{M, T}(Tuple(T.(conf[:u_bg])))
    u_box = State{M, T}(Tuple(T.(conf[:u_box])))
    mins = Space{D, T}(Tuple(T.(conf[:mins])))
    maxs = Space{D, T}(Tuple(T.(conf[:maxs])))
    
    return Box(u_bg, u_box, mins, maxs)
end