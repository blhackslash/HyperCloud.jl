struct Riemann{D, M} <: InitialCondition
    uL::State{M}
    uR::State{M}
    p0::Space{D}
    n::Space{D}
end

(ic::Riemann)(pos::Space{D}) where {D} = dot(pos - ic.p0, ic.n) < 0 ? ic.uL : ic.uR

function build_ic(::Val{:riemann}, p, D::Int, ::Type{T}) where {T}
    uL = param2uvec(p[1])
    uR = param2uvec(p[2])
    p0 = param2xvec(p[3])
    # If normal vector not provided, default to pointing in +X direction
    n = length(p) > 3 ? normalize(param2xvec(p[4])) : SVector{length(p0), Float64}(ntuple(i -> i==1 ? 1.0 : 0.0, length(p0)))
    
    return Riemann(uL, uR, p0, n)
end