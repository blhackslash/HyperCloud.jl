struct Box{D, M} <: InitialCondition
    u_bg::State{M}
    u_box::State{M}
    mins::Space{D}
    maxs::Space{D}
end

(ic::Box)(pos::Space{D}) where {D} = all(ic.mins .<= pos .<= ic.maxs) ? ic.u_box : ic.u_bg

function build_ic(::Val{:box}, p, D::Int, ::Type{T}) where {T}
    if length(p) == 4
        # 1D/Multi-D unified format
        return Box(param2uvec(p[1]), param2uvec(p[2]), param2xvec(p[3]), param2xvec(p[4]))
    elseif length(p) == 6
        # Backwards compatibility for old 2D format
        mins = param2xvec((p[3], p[5]))
        maxs = param2xvec((p[4], p[6]))
        return Box(param2uvec(p[1]), param2uvec(p[2]), mins, maxs)
    else
        error("Invalid number of parameters for Box")
    end
end