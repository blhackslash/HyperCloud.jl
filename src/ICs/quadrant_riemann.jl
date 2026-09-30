struct QuadrantRiemann{D, M, N_states} <: InitialCondition
    u_states::NTuple{N_states, State{M}} 
    p0::Space{D}
end

function (ic::QuadrantRiemann{D, M, N})(pos::Space{D}) where {D, M, N}
    # Binary encoding: Left/Bottom adds 0, Right/Top adds 2^(d-1)
    idx = 1
    for d in 1:D
        if pos[d] >= ic.p0[d]
            idx += 2^(d-1)
        end
    end
    return ic.u_states[idx]
end

function build_ic(::Val{:q_riemann}, p, D::Int, ::Type{T}) where {T}
    states = ntuple(i -> param2uvec(p[1][i]), length(p[1]))
    p0 = param2xvec(p[2])
    return QuadrantRiemann(states, p0)
end