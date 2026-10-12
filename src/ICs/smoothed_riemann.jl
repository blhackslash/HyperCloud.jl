export SRiemann

"""
    SRiemann{D, M, T} <: InitialCondition

Smoothed, continuous approximation of a standard Riemann problem that avoids infinite initial gradients at the interface.

# Fields
- `uL::State{M, T}`: The asymptotic physical state on the "left".
- `uR::State{M, T}`: The asymptotic physical state on the "right".
- `p0::Space{D, T}`: A point located on the central interface.
- `n::Space{D, T}`: The normal vector pointing toward the right state.
- `width::T`: Controls the steepness of the transition.

Uses an arctangent to blend smoothly between `uL` and `uR` over the given `width`.
All fields are converted to the precision `T` of `uL`.
"""
struct SRiemann{D, M, T} <: InitialCondition
    uL::State{M, T}
    uR::State{M, T}
    p0::Space{D, T}
    n::Space{D, T}
    width::T
end

# Mixed-precision convenience constructor (the default one requires matching T everywhere)
SRiemann(uL::State{M, T}, uR::State{M}, p0::Space{D}, n::Space{D}, width::Real) where {D, M, T} =
    SRiemann{D, M, T}(uL, uR, p0, n, width)

function (ic::SRiemann{D, M, T})(pos::Space{D}) where {D, M, T}
    dist = dot(pos - ic.p0, ic.n)
    return (ic.uL + ic.uR) / 2 - (ic.uL - ic.uR) * (atan(dist / ic.width) / T(π))
end

function build_ic(::Val{:s_riemann}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int

    uL = State{M, T}(Tuple(T.(conf[:uL])))
    uR = State{M, T}(Tuple(T.(conf[:uR])))
    p0 = Space{D, T}(Tuple(T.(conf[:p0])))
    width = T(conf[:width])

    n_tup = get(conf, :n, ntuple(i -> i == 1 ? 1.0 : 0.0, D))
    n = normalize(Space{D, T}(Tuple(T.(n_tup))))

    return SRiemann(uL, uR, p0, n, width)
end
