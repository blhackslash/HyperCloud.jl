export BurgersEquation

"""
    BurgersEquation{D, T} <: HyperbolicPDE{D, 1, T}

Implements the multi-dimensional inviscid Burgers' equation. 
The system mathematically enforces a strictly scalar state vector (`M = 1`) across any spatial dimension `D`.
"""
struct BurgersEquation{D, T} <: HyperbolicPDE{D, 1, T} end

@inline prim2cons(::BurgersEquation, u::State) = u
@inline cons2prim(::BurgersEquation, w::State) = w

@inline function flux(eq::BurgersEquation{D, T}, u::State{1, T}) where {D, T}
    return Flux{D, 1, T}(ntuple(_ -> State{1, T}(T(0.5) * u[1]^2), Val(D)))
end

@inline function max_eigenvalue(eq::BurgersEquation, u::State{1, T}, d::Int) where {T}
    return abs(u[1])
end

# Returns the localized Jacobian (u) as a strict 1x1 SMatrix
@inline velocity(eq::BurgersEquation, u::State{1, T}, d::Int) where {T} = SMatrix{1, 1, T, 1}(u[1])

# Moving grids: transport velocity (u in every direction) and exact relative wave speed
@inline fluid_velocity(eq::BurgersEquation{D, T}, u::State{1, T}) where {D, T} = Space{D, T}(ntuple(_ -> u[1], Val(D)))
@inline max_eigenvalue(eq::BurgersEquation, u::State{1, T}, d::Int, shift::Real) where {T} = abs(u[1] - T(shift))

# =========================================================================
# BURGERS ANALYTIC CLOSURES
# =========================================================================

# Fallback for ICs without an exact solution: returns the initial condition, which is only
# correct at t = 0. Warn once at construction so a wrong reference is never silent.
function analytic_closure(eq::BurgersEquation{D, T}, ic::InitialCondition, geom::GeometricDomain) where {D, T}
    @warn "No exact Burgers solution for $(nameof(typeof(ic))) in $(D)D; the reference returns the initial condition and is only valid at t = 0."
    return (st::SVector) -> ic(SVector{D, T}(ntuple(d -> T(st[d]), Val(D))))
end

# Multi-Dimensional Gauss (Newton-Raphson Characteristic Trace)
function analytic_closure(eq::BurgersEquation{D, T}, ic::Gauss, geom::GeometricDomain) where {D, T}
    tol = T(1e-10)
    max_iter = 100
    
    return function exact_burgers_gauss(st::SVector)
        t = T(st[end])
        pos = SVector{D, T}(ntuple(d -> T(st[d]), Val(D)))
        
        if t <= T(1e-12); return ic(pos); end
        
        a = ic.a[1]
        b = ic.b
        w = ic.width
        w2 = w^2
        
        u_curr = ic(pos)[1]
        
        for _ in 1:max_iter
            x0 = pos .- u_curr .* t
            dist2 = sum(abs2, x0 .- b)
            exp_term = exp(-dist2 / w2)
            
            F = u_curr - a * exp_term
            sum_diff = sum(x0 .- b)
            dFdu = one(T) - (a * exp_term) * (T(2.0) * t * sum_diff / w2)
            
            u_next = u_curr - F / dFdu
            
            if abs(u_next - u_curr) < tol
                return SVector{1, T}(u_next)
            end
            u_curr = u_next
        end
        
        @warn "Newton-Raphson did not converge for ND Burgers Gauss at pos=$pos, t=$t."
        return SVector{1, T}(u_curr)
    end
end

# Planar Riemann problem in any dimension D.
#
# With ξ = (x - p0)⋅n and u = u(ξ, t), the equation u_t + Σ_d ∂_d(u²/2) = 0 becomes
# u_t + a ∂_ξ(u²/2) = 0 with a = Σ_d n_d. The wave type therefore depends on the
# characteristic speeds a·uL and a·uR, not on uL vs uR alone: flipping the sign of `n`
# (and swapping uL/uR) describes the same data and must give the same solution.
# `n` does not need to be normalized.
#
# For SRiemann this returns the sharp-interface Riemann solution (the t → 0 limit of the
# smoothed data is ignored), which is a good reference once the fan/shock is wider than `width`.
function analytic_closure(eq::BurgersEquation{D, T}, ic::Union{Riemann, SRiemann}, geom::GeometricDomain) where {D, T}
    p0 = SVector{D, T}(ic.p0)
    n  = SVector{D, T}(ic.n)
    uL, uR = T(ic.uL[1]), T(ic.uR[1])

    a = sum(n)
    aL, aR = a * uL, a * uR      # characteristic speeds along ξ
    s = a * (uL + uR) / 2        # shock speed along ξ

    return function exact_burgers_riemann(st::SVector)
        t = T(st[end])
        pos = SVector{D, T}(ntuple(d -> T(st[d]), Val(D)))

        if t <= T(1e-12); return SVector{1, T}(ic(pos)); end

        ξ = dot(pos - p0, n)
        if aL > aR                         # converging characteristics ⇒ shock
            return SVector{1, T}(ξ < s * t ? uL : uR)
        elseif ξ < aL * t                  # left of the fan
            return SVector{1, T}(uL)
        elseif ξ > aR * t                  # right of the fan
            return SVector{1, T}(uR)
        elseif abs(a) < T(1e-14)           # stationary interface (a = 0): same convention as the IC
            return SVector{1, T}(uR)
        else                               # inside the fan
            return SVector{1, T}(ξ / (a * t))
        end
    end
end

# 1D Box
function analytic_closure(eq::BurgersEquation{1, T}, ic::Box, geom::GeometricDomain) where {T}
    xs, xe = ic.mins[1], ic.maxs[1]
    ub, ug = ic.u_box[1], ic.u_bg[1]
    
    return function exact_burgers_box(st::SVector)
        t = T(st[end])
        x = T(st[1])
        
        if t <= T(1e-12); return ic(SVector{1, T}(x)); end
        if abs(ub - ug) < T(1e-12); return ic.u_bg; end

        if ub > ug # Top-hat case
            t_int = T(2.0) * (xe - xs) / (ub - ug)
            if t < t_int
                s_shock = T(0.5) * (ub + ug)
                if x < xs + ug * t; return ic.u_bg
                elseif x < xs + ub * t; return SVector{1, T}((x - xs) / t)
                elseif x < xe + s_shock * t; return ic.u_box
                else return ic.u_bg
                end
            else
                C = sqrt(T(2.0) * (xe - xs) * (ub - ug))
                x_shock = xs + ug * t + C * sqrt(t)
                if x < xs + ug * t; return ic.u_bg
                elseif x < x_shock; return SVector{1, T}((x - xs) / t)
                else return ic.u_bg
                end
            end
        else # Well case
            t_int = T(2.0) * (xe - xs) / (ug - ub)
            if t < t_int
                s_shock = T(0.5) * (ug + ub)
                if x < xs + s_shock * t; return ic.u_bg
                elseif x < xe + ub * t; return ic.u_box
                elseif x < xe + ug * t; return SVector{1, T}((x - xe) / t)
                else return ic.u_bg
                end
            else
                C = sqrt(T(2.0) * (xe - xs) * (ug - ub))
                x_shock = xe + ug * t - C * sqrt(t)
                if x < x_shock; return ic.u_bg
                elseif x < xe + ug * t; return SVector{1, T}((x - xe) / t)
                else return ic.u_bg
                end
            end
        end
    end
end
