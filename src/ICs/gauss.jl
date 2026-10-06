export Gauss

"""
    Gauss{D, M} <: InitialCondition

Evaluates a smooth, multi-dimensional Gaussian pulse.

# Fields
- `a::State{M}`: The amplitude state vector.
- `b::Space{D}`: The spatial center coordinate of the pulse.
- `width::Float64`: The parameter controlling the spatial spread of the bell curve.

The physical state is computed mathematically as `a * exp(-sum(abs2, pos - b) / width^2)`.
"""
struct Gauss{D, M} <: InitialCondition
    a::State{M}
    b::Space{D}
    width::Float64
end

(ic::Gauss)(pos::Space{D}) where {D} = ic.a * exp(-sum(abs2, pos - ic.b) / ic.width^2)

function build_ic(::Val{:gauss}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int
    
    a = State{M, T}(Tuple(T.(conf[:a])))
    b = Space{D, T}(Tuple(T.(conf[:b])))
    width = T(conf[:width])
    
    return Gauss(a, b, width)
end