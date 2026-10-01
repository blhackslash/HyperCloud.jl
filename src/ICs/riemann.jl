struct Riemann{D, M} <: InitialCondition
    uL::State{M}
    uR::State{M}
    p0::Space{D}
    n::Space{D}
end

(ic::Riemann)(pos::Space{D}) where {D} = dot(pos - ic.p0, ic.n) < 0 ? ic.uL : ic.uR

function build_ic(::Val{:riemann}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:Type]::DataType, ctx[:D]::Int, ctx[:M]::Int
    
    uL = State{M, T}(Tuple(T.(conf[:uL])))
    uR = State{M, T}(Tuple(T.(conf[:uR])))
    p0 = Space{D, T}(Tuple(T.(conf[:p0])))
    
    # Default normal points in +X direction if not provided
    n_tup = get(conf, :n, ntuple(i -> i == 1 ? 1.0 : 0.0, D))
    n = normalize(Space{D, T}(Tuple(T.(n_tup))))
    
    return Riemann(uL, uR, p0, n)
end