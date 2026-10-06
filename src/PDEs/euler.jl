export EulerEquation

"""
    EulerEquation{D, M, T, R} <: HyperbolicPDE{D, M, T}

Implements the compressible Euler equations for ideal gases. 

# Fields
- `gamma::T`: The specific heat ratio of the gas.
- `rep::R`: The equation representation (e.g., `Conservative` or `NCRepresentation`).

# Constructors
- `EulerEquation(::Val{D}, ::Type{T}, gamma::T = T(1.4), rep::R = Conservative())`: Automatically configures the system for spatial dimension `D` and sets the state vector length to `M = D + 2`.
"""
struct EulerEquation{D, M, T, R} <: HyperbolicPDE{D, M, T}
    gamma::T
    rep::R
end

function EulerEquation(::Val{D}, ::Type{T}, gamma::T = T(1.4), rep::R = Conservative()) where {D, T, R <: EquationRepresentation}
    return EulerEquation{D, D + 2, T, R}(gamma, rep)
end

@inline function prim2cons(eq::EulerEquation{D, M, T}, V::State{M, T}) where {D, M, T}
    rho = V[1]
    u = Space{D, T}(ntuple(d -> V[1+d], Val(D))) 
    p = V[M]
    
    m = rho .* u
    E = p / (eq.gamma - T(1)) + T(0.5) * rho * sum(abs2, u) 
    
    return State{M, T}(rho, m..., E)
end

@inline function cons2prim(eq::EulerEquation{D, M, T}, U::State{M, T}) where {D, M, T}
    rho = max(U[1], T(1e-7))
    m = Space{D, T}(ntuple(d -> U[1+d], Val(D)))
    E = U[M]
    
    u = m ./ rho
    p = (eq.gamma - T(1)) * (E - T(0.5) * sum(abs2, m) / rho)
    
    return State{M, T}(rho, u..., max(p, T(1e-7)))
end

@inline function flux(eq::EulerEquation{D, M, T, <:Conservative}, U::State{M, T}) where {D, M, T}
    V = cons2prim(eq, U)
    rho = V[1]
    u = Space{D, T}(ntuple(d -> V[1+d], Val(D)))
    p = V[M]
    E = U[M]
    
    return Flux{D, M, T}(ntuple(Val(D)) do d
        ud = u[d]
        mass_flux = rho * ud
        mom_flux = ntuple(i -> rho * u[i] * ud + (i == d ? p : zero(T)), Val(D))
        energy_flux = ud * (E + p)
        
        State{M, T}(mass_flux, mom_flux..., energy_flux)
    end)
end

@inline function max_eigenvalue(eq::EulerEquation{D, M, T, <:Conservative}, U::State{M, T}, d::Int) where {D, M, T}
    rho = max(U[1], T(1e-9))
    m_d = U[1+d]
    E = U[M]
    
    m_sq = sum(abs2, ntuple(i -> U[1+i], Val(D)))
    p = max((eq.gamma - one(T)) * (E - T(0.5) * m_sq / rho), T(1e-9))
    c = sqrt(eq.gamma * p / rho)
    
    return abs(m_d / rho) + c
end

@inline function max_eigenvalue(eq::EulerEquation{D, M, T, <:NCRepresentation}, V::State{M, T}, d::Int) where {D, M, T}
    rho = max(V[1], T(1e-9))
    u_d = V[1+d]
    p = max(V[M], T(1e-9))
    
    c = sqrt(eq.gamma * p / rho)
    
    return abs(u_d) + c
end

# =========================================================================
# EULER ANALYTIC CLOSURES
# =========================================================================

function analytic_closure(eq::EulerEquation{1, M, T}, ic::Riemann, geom::GeometricDomain) where {M, T}
    gamma = eq.gamma
    
    rho_L, m_L, E_L = ic.uL
    rho_R, m_R, E_R = ic.uR
    x0 = ic.p0[1]
    
    u_L = m_L / rho_L
    p_L = (gamma - one(T)) * (E_L - T(0.5) * rho_L * u_L^2)
    u_R = m_R / rho_R
    p_R = (gamma - one(T)) * (E_R - T(0.5) * rho_R * u_R^2)
    
    c_L = sqrt(gamma * p_L / rho_L)
    c_R = sqrt(gamma * p_R / rho_R)
    
    # Pre-calculate exact p_star and u_star
    function pressure_func(p_star_guess::T)
        f_L = p_star_guess > p_L ? (p_star_guess - p_L) * sqrt((T(2.0) / ((gamma + one(T)) * rho_L)) / (p_star_guess + p_L * (gamma - one(T)) / (gamma + one(T)))) : (T(2.0) * c_L / (gamma - one(T))) * ((p_star_guess / p_L)^((gamma - one(T)) / (T(2.0) * gamma)) - one(T))
        f_R = p_star_guess > p_R ? (p_star_guess - p_R) * sqrt((T(2.0) / ((gamma + one(T)) * rho_R)) / (p_star_guess + p_R * (gamma - one(T)) / (gamma + one(T)))) : (T(2.0) * c_R / (gamma - one(T))) * ((p_star_guess / p_R)^((gamma - one(T)) / (T(2.0) * gamma)) - one(T))
        return f_L + f_R + (u_R - u_L)
    end

    p_star = T(0.5) * (p_L + p_R)
    for _ in 1:100
        f_p = pressure_func(p_star)
        if abs(f_p) < T(1e-9); break; end
        dfdp = (pressure_func(p_star * T(1.001)) - f_p) / (p_star * T(0.001))
        p_star = max(T(1e-9), p_star - f_p / (dfdp + T(1e-9)))
    end

    f_L_final = p_star > p_L ? (p_star - p_L) * sqrt((T(2.0) / ((gamma + one(T)) * rho_L)) / (p_star + p_L * (gamma - one(T)) / (gamma + one(T)))) : (T(2.0) * c_L / (gamma - one(T))) * ((p_star / p_L)^((gamma - one(T)) / (T(2.0) * gamma)) - one(T))
    u_star = u_L - f_L_final

    return function exact_euler_riemann(st::SVector)
        t = T(st[end])
        x = T(st[1])
        
        if t <= T(1e-9); return ic(SVector{1, T}(x)); end
        
        s_query = (x - x0) / t
        rho_final, u_final, p_final = zero(T), zero(T), zero(T)

        if s_query <= u_star # Left of contact
            if p_star > p_L # Left Shock
                S_L = u_L - c_L * sqrt((gamma + one(T)) / (T(2.0) * gamma) * (p_star / p_L) + (gamma - one(T)) / (T(2.0) * gamma))
                rho_star_L = rho_L * ((p_star / p_L) + (gamma - one(T)) / (gamma + one(T))) / (one(T) + (p_star / p_L) * (gamma - one(T)) / (gamma + one(T)))
                rho_final, u_final, p_final = s_query <= S_L ? (rho_L, u_L, p_L) : (rho_star_L, u_star, p_star)
            else # Left Rarefaction
                S_head_L = u_L - c_L
                S_tail_L = u_star - c_L * (p_star / p_L)^((gamma - one(T)) / (T(2.0) * gamma))
                if s_query <= S_head_L; rho_final, u_final, p_final = rho_L, u_L, p_L
                elseif s_query >= S_tail_L; rho_final, u_final, p_final = rho_L * (p_star / p_L)^(one(T) / gamma), u_star, p_star
                else
                    u_final = (T(2.0) / (gamma + one(T))) * (c_L + (gamma - one(T)) / T(2.0) * u_L + s_query)
                    c_final = c_L - (gamma - one(T)) / T(2.0) * (u_final - u_L)
                    rho_final = rho_L * (c_final / c_L)^(T(2.0) / (gamma - one(T)))
                    p_final = p_L * (rho_final / rho_L)^gamma
                end
            end
        else # Right of contact
            if p_star > p_R # Right Shock
                S_R = u_R + c_R * sqrt((gamma + one(T)) / (T(2.0) * gamma) * (p_star / p_R) + (gamma - one(T)) / (T(2.0) * gamma))
                rho_star_R = rho_R * ((p_star / p_R) + (gamma - one(T)) / (gamma + one(T))) / (one(T) + (p_star / p_R) * (gamma - one(T)) / (gamma + one(T)))
                rho_final, u_final, p_final = s_query >= S_R ? (rho_R, u_R, p_R) : (rho_star_R, u_star, p_star)
            else # Right Rarefaction
                S_head_R = u_R + c_R
                S_tail_R = u_star + c_R * (p_star / p_R)^((gamma - one(T)) / (T(2.0) * gamma))
                if s_query >= S_head_R; rho_final, u_final, p_final = rho_R, u_R, p_R
                elseif s_query <= S_tail_R; rho_final, u_final, p_final = rho_R * (p_star / p_R)^(one(T) / gamma), u_star, p_star
                else
                    u_final = (T(2.0) / (gamma + one(T))) * (-c_R + (gamma - one(T)) / T(2.0) * u_R + s_query)
                    c_final = c_R + (gamma - one(T)) / T(2.0) * (u_final - u_R)
                    rho_final = rho_R * (c_final / c_R)^(T(2.0) / (gamma - one(T)))
                    p_final = p_R * (rho_final / rho_R)^gamma
                end
            end
        end

        return SVector{3, T}(rho_final, rho_final * u_final, p_final / (gamma - one(T)) + T(0.5) * rho_final * u_final^2)
    end
end