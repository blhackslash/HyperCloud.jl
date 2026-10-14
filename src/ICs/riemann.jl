export Riemann

"""
    Riemann{D, M, T} <: InitialCondition

Defines a standard discontinuous Riemann problem (e.g., a shock tube or half-space) split by a hyperplane.

# Fields
- `uL::State{M, T}`: The physical state vector on the "left" side of the hyperplane.
- `uR::State{M, T}`: The physical state vector on the "right" side of the hyperplane.
- `p0::Space{D, T}`: An origin point located exactly on the separating hyperplane.
- `n::Space{D, T}`: The normal vector of the hyperplane pointing toward the right state.

Returns `uL` if the dot product of the relative position and the normal is negative, and `uR` otherwise.
"""
struct Riemann{D, M, T} <: InitialCondition
    uL::State{M, T}
    uR::State{M, T}
    p0::Space{D, T}
    n::Space{D, T}
end

function Riemann(uL::State{M}, uR::State{M}, p0::Space{D}, n::Space{D}) where {D, M}
    T = promote_type(eltype(uL), eltype(uR), eltype(p0), eltype(n))
    return Riemann{D, M, T}(uL, uR, p0, n)
end

(ic::Riemann)(pos::Space{D}) where {D} = dot(pos - ic.p0, ic.n) < 0 ? ic.uL : ic.uR

function build_ic(::Val{:riemann}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int
    
    uL = State{M, T}(Tuple(T.(conf[:uL])))
    uR = State{M, T}(Tuple(T.(conf[:uR])))
    p0 = Space{D, T}(Tuple(T.(conf[:p0])))
    
    # Default normal points in +X direction if not provided
    n_tup = get(conf, :n, ntuple(i -> i == 1 ? 1.0 : 0.0, D))
    n = normalize(Space{D, T}(Tuple(T.(n_tup))))
    
    return Riemann(uL, uR, p0, n)
end