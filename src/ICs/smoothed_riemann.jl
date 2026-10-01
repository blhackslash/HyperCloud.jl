struct SRiemann{D, M} <: SmoothInitialCondition
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
    T, D, M = ctx[:Type]::DataType, ctx[:D]::Int, ctx[:M]::Int
    
    uL = State{M, T}(Tuple(T.(conf[:uL])))
    uR = State{M, T}(Tuple(T.(conf[:uR])))
    p0 = Space{D, T}(Tuple(T.(conf[:p0])))
    width = T(conf[:width])
    
    n_tup = get(conf, :n, ntuple(i -> i == 1 ? 1.0 : 0.0, D))
    n = normalize(Space{D, T}(Tuple(T.(n_tup))))
    
    return SRiemann(uL, uR, p0, n, width)
end