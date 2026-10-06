export Box

"""
    Box{D, M} <: InitialCondition

Defines a discontinuous rectangular sub-region.

# Fields
- `u_bg::State{M}`: The background state applied outside the box.
- `u_box::State{M}`: The interior state applied within the box.
- `mins::Space{D}`: The lower spatial bounds of the rectangular region.
- `maxs::Space{D}`: The upper spatial bounds of the rectangular region.

Evaluates to `u_box` if the queried position falls strictly within `mins` and `maxs`, and evaluates to `u_bg` otherwise.
"""
struct Box{D, M} <: InitialCondition
    u_bg::State{M}
    u_box::State{M}
    mins::Space{D}
    maxs::Space{D}
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