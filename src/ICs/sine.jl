struct Sine{D, M} <: SmoothInitialCondition
    a::State{M}
    period::Space{D}
    c_offset::State{M}
end

(ic::Sine)(pos::Space{D}) where {D} = ic.a * sin(2.0 * pi * sum(pos ./ ic.period)) + ic.c_offset

function build_ic(::Val{:sine}, p, D::Int, ::Type{T}) where {T}
    return Sine(param2uvec(p[1]), param2xvec(p[2]), param2uvec(p[3]))
end