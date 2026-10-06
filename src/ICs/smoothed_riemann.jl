export SRiemann

"""
    SRiemann{D, M} <: InitialCondition

Provides a smoothed, continuous approximation of a standard Riemann problem to prevent infinite initial gradients at the interface[cite: 32].

# Fields
- `uL::State{M}`: The asymptotic physical state on the "left"[cite: 32].
- `uR::State{M}`: The asymptotic physical state on the "right"[cite: 32].
- `p0::Space{D}`: A point located on the central interface[cite: 32].
- `n::Space{D}`: The normal vector pointing toward the right state[cite: 32].
- `width::Float64`: The width parameter controlling the steepness of the continuous transition[cite: 32].

Utilizes an arctangent function to transition smoothly between `uL` and `uR` over the defined `width`[cite: 32].
"""
struct SRiemann{D, M} <: InitialCondition
    uL::State{M}
    uR::State{M}
    p0::Space{D}
    n::Space{D}
    width::Float64
end

function (ic::SRiemann)(pos::Space{D}) where {D}
    dist = dot(pos - ic.p0, ic.n)
    return @. 0.5 * (ic.uL + ic.uR) - (ic.uL - ic.uR) / pi * atan(dist / ic.width)
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