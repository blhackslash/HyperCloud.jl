struct Gauss{D, M} <: SmoothInitialCondition
    a::State{M}
    b::Space{D}
    width::Float64
end

(ic::Gauss)(pos::Space{D}) where {D} = ic.a * exp(-sum(abs2, pos - ic.b) / ic.width^2)

function build_ic(::Val{:gauss}, p, D::Int, ::Type{T}) where {T}
    return Gauss(param2uvec(p[1]), param2xvec(p[2]), Float64(p[3]))
end