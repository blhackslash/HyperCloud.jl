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

function build_ic(::Val{:s_riemann}, p, D::Int, ::Type{T}) where {T}
    uL = param2uvec(p[1])
    uR = param2uvec(p[2])
    p0 = param2xvec(p[3])
    
    if length(p) == 4
        n = SVector{length(p0), Float64}(ntuple(i -> i==1 ? 1.0 : 0.0, length(p0)))
        return SRiemann(uL, uR, p0, n, Float64(p[4]))
    else
        n = normalize(param2xvec(p[4]))
        return SRiemann(uL, uR, p0, n, Float64(p[5]))
    end
end