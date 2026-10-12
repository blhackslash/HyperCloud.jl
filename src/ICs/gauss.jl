export Gauss

"""
    Gauss{D, M, T} <: InitialCondition

Evaluates a smooth, multi-dimensional Gaussian pulse.

# Fields
- `a::State{M, T}`: The amplitude state vector.
- `b::Space{D, T}`: The spatial center coordinate of the pulse.
- `width::T`: The parameter controlling the spatial spread of the bell curve.

The physical state is computed as `a * exp(-sum(abs2, pos - b) / width^2)`.
The center and width are converted to the precision `T` of the amplitude, so a
`Float32` amplitude gives a pure `Float32` pulse even if `b` or `width` are `Float64`.
"""
struct Gauss{D, M, T} <: InitialCondition
    a::State{M, T}
    b::Space{D, T}
    width::T
end

# Mixed-precision convenience constructor (the default one requires matching T everywhere)
Gauss(a::State{M, T}, b::Space{D}, width::Real) where {D, M, T} = Gauss{D, M, T}(a, b, width)

(ic::Gauss)(pos::Space{D}) where {D} = ic.a * exp(-sum(abs2, pos - ic.b) / ic.width^2)

function build_ic(::Val{:gauss}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int

    a = State{M, T}(Tuple(T.(conf[:a])))
    b = Space{D, T}(Tuple(T.(conf[:b])))
    width = T(conf[:width])

    return Gauss(a, b, width)
end
