struct QRiemann{D, M, N_states} <: InitialCondition
    u_states::NTuple{N_states, State{M}} 
    p0::Space{D}
end

function (ic::QRiemann{D, M, N})(pos::Space{D}) where {D, M, N}
    # Binary encoding: Left/Bottom adds 0, Right/Top adds 2^(d-1)
    idx = 1
    for d in 1:D
        if pos[d] >= ic.p0[d]
            idx += 2^(d-1)
        end
    end
    return ic.u_states[idx]
end

function build_ic(::Val{:q_riemann}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:Type]::DataType, ctx[:D]::Int, ctx[:M]::Int
    
    # Process the tuple of states explicitly via the extracted dictionary list
    raw_states = conf[:states]::Tuple
    states = ntuple(i -> State{M, T}(Tuple(T.(raw_states[i]))), length(raw_states))
    
    p0 = Space{D, T}(Tuple(T.(conf[:p0])))
    
    return QRiemann(states, p0)
end