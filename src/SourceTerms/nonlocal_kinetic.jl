export NonLocalRelaxationSourceTerm

mutable struct NonLocalRelaxationSourceTerm{D, NM, NK, T, MEQ <: HyperbolicPDE} <: AbstractImplicitSourceTerm
    km::Kin2Macro{NM, NK}
    macro_eq::MEQ
    inv_epsilon::T
    coefficients::State{NM, T}
    scaled_inv_speeds::SVector{NK, Space{D, T}}
    t_potential::Matrix{T}
end

function NonLocalRelaxationSourceTerm(
    km::Kin2Macro{NM, NK}, 
    eps::T, 
    coeffs::State{NM, T}, 
    macro_eq::HyperbolicPDE{D, NM, T},
    eq_kin::HyperbolicPDE{D, NK, T}, 
    interior_factor::T = T(D)
) where {D, NM, NK, T}
    
    scaled_inv_speeds = ntuple(Val(NK)) do k
        Space{D, T}(ntuple(Val(D)) do d
            v = kinetic_wave_speed(eq_kin, d, k)
            abs(v) > T(1e-14) ? interior_factor / v : zero(T)
        end)
    end
    
    t_potential = Matrix{T}(undef, 0, NM)
    return NonLocalRelaxationSourceTerm(
        km, macro_eq, one(T) / eps, coeffs, SVector{NK, Space{D, T}}(scaled_inv_speeds), t_potential
    )
end

function ensure_buffer_size!(st::NonLocalRelaxationSourceTerm{D, NM, NK, T}, N_particles::Int) where {D, NM, NK, T}
    if size(st.t_potential, 1) != N_particles
        st.t_potential = Matrix{T}(undef, N_particles, NM)
    end
end

function pre_solve_update!(st::NonLocalRelaxationSourceTerm{D, NM, NK, T}, stage_data::AbstractVector{State{NK, T}}, pg::ParticleGrid, t::Real) where {D, NM, NK, T}
    N_particles = pg.meta.N
    ensure_buffer_size!(st, N_particles)
    
    km = st.km
    macro_eq = st.macro_eq
    
    # 1. Compute jump potentials iteratively along the array
    Threads.@threads for i in 2:N_particles
        v_L = stage_data[i-1]
        v_R = stage_data[i]    
        
        u_L = km(v_L)
        u_R = km(v_R)
        
        jump = path_integral(macro_eq, u_L, u_R)
        
        for m in 1:NM
            st.t_potential[i, m] = jump[m]
        end
    end
    
    # Base condition for particle 1
    for m in 1:NM; st.t_potential[1, m] = zero(T); end
    
    # 2. Cumulative summation to build the topological field
    for i in 2:N_particles
        for m in 1:NM
            st.t_potential[i, m] += st.t_potential[i-1, m]
        end
    end
    return nothing
end

@inline function evaluate_source(st::NonLocalRelaxationSourceTerm{D, NM, NK, T}, V_kin::State{NK, T}, p_idx::Int, pg::ParticleGrid, t::Real) where {D, NM, NK, T}
    u_macro = st.km(V_kin)
    
    return State{NK, T}(ntuple(Val(NK)) do k
        m_idx = st.km(k)
        T_val = st.t_potential[p_idx, m_idx]
        
        T_dot_inv_lambda = T_val * st.scaled_inv_speeds[k][1]
        Mk_val = st.coefficients[m_idx] * (u_macro[m_idx] + T_dot_inv_lambda)
        
        (Mk_val - V_kin[k]) * st.inv_epsilon
    end)
end

# --- NON-LOCAL Relaxation Initialization ---
function set_initial_conditions!(
    pg::ParticleGrid{D, NK},
    st::NonLocalRelaxationSourceTerm{D, NM, NK},
    IC::InitialCondition,
    eq_macro::HyperbolicPDE{D}
) where {D, NM, NK}
    
    N_particles = pg.meta.N
    
    # 1. Initialize grid to LOCAL equilibrium (V_k = c_m * U_m)
    for p_idx in 1:N_particles
        u_val = IC(pg.core.positions[p_idx]) 
        pg.rhos[p_idx] = State{NK}(ntuple(Val(NK)) do k
            st.coefficients[st.km(k)] * u_val[st.km(k)]
        end)
    end

    # 2. Compute the true initial potential T_0 using current grid state
    update_nonlocal_potential!(st, pg.rhos, pg, eq_macro)

    # 3. Re-initialize kinetic grids to the NON-LOCAL equilibrium: V_0 = M(U_0, T_0)
    for p_idx in 1:N_particles
        u_val = IC(pg.core.positions[p_idx])
        
        pg.rhos[p_idx] = State{NK}(ntuple(Val(NK)) do k
            m_idx = st.km(k)
            T_val = st.T_potential[p_idx, m_idx]
            T_dot_inv_lambda = T_val * st.scaled_inv_speeds[k][1]
            
            return st.coefficients[m_idx] * (u_val[m_idx] + T_dot_inv_lambda)
        end)
    end
    
    @info "Initialized Non-Local Equilibrium (Max Potential: $(maximum(abs.(st.T_potential))))"
    return nothing
end
