struct BurgersEquation{D, T, R} <: HyperbolicPDE{D, 1, T, R} 
    rep::R
end

BurgersEquation(::Val{D}, ::Type{T}, rep::R = Conservative()) where {D, T, R <: EquationRepresentation} = BurgersEquation{D, T, R}(rep)

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

function build_equation(::Val{:burgers}, pde_conf::Dict, context::Dict)
    T = context[:Type]::DataType
    D = pde_conf[:D]::Int
    rep = parse_representation(pde_conf)
    
    return BurgersEquation(Val(D), T, rep)
end

# Generic fallback for t=0 or unhandled ICs
analytic_closure(eq::BurgersEquation{D}, ic::InitialCondition, params::ParamDict) where {D} = (st::SVector) -> ic(SVector{D, Float64}(ntuple(d -> st[d], Val(D))))

# =========================================================================
# BURGERS ANALYTIC CLOSURES
# =========================================================================

# Generic fallback for t=0 or unhandled ICs
analytic_closure(eq::BurgersEquation{D, T}, ic::InitialCondition, geom::GeometricDomain) where {D, T} = 
    (st::SVector) -> ic(SVector{D, T}(ntuple(d -> T(st[d]), Val(D))))

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

# 1D Sine
function analytic_closure(eq::BurgersEquation{1, T}, ic::Sine, geom::GeometricDomain) where {T}
    tol = T(1e-10)
    max_iter = 100
    return function exact_burgers_sine(st::SVector)
        t = T(st[end])
        x = T(st[1])
        pos = SVector{1, T}(x)
        
        if t <= T(1e-12); return ic(pos); end
        
        u_curr = ic(pos)[1]
        
        for _ in 1:max_iter
            u_next = ic.a[1] * sin(T(2.0 * pi) * (x - u_curr * t) / ic.period[1]) + ic.c_offset[1]
            if abs(u_next - u_curr) < tol; return SVector{1, T}(u_next); end
            u_curr = u_next
        end
        return SVector{1, T}(u_curr)
    end
end

# 1D Riemann
function analytic_closure(eq::BurgersEquation{1, T}, ic::Union{Riemann, SRiemann}, geom::GeometricDomain) where {T}
    x0 = ic.p0[1]
    uL = ic.uL[1]
    uR = ic.uR[1]
    s = T(0.5) * (uL + uR)

    return function exact_burgers_riemann(st::SVector)
        t = T(st[end])
        x = T(st[1])
        
        if t <= T(1e-12); return ic(SVector{1, T}(x)); end
        
        if uL > uR # Shock
            return x < x0 + s * t ? ic.uL : ic.uR
        else # Rarefaction
            if x < x0 + uL * t; return ic.uL
            elseif x > x0 + uR * t; return ic.uR
            else return SVector{1, T}((x - x0) / t)
            end
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

# 2D Riemann
function analytic_closure(eq::BurgersEquation{2, T}, ic::Riemann, geom::GeometricDomain) where {T}
    n_sum = ic.n[1] + ic.n[2]
    uL, uR = ic.uL[1], ic.uR[1]
    s = T(0.5) * (uL + uR) * n_sum

    return function exact_burgers2d_riemann(st::SVector)
        t = T(st[end])
        pos = SVector{2, T}(T(st[1]), T(st[2]))
        
        if t <= T(1e-12); return ic(pos); end
        
        d = dot(pos - ic.p0, ic.n)
        if uL > uR # Shock
            return d < s * t ? ic.uL : ic.uR
        else # Rarefaction
            if d < uL * n_sum * t; return ic.uL
            elseif d > uR * n_sum * t; return ic.uR
            else 
                if abs(t * n_sum) < T(1e-14); return SVector{1, T}(T(0.5) * (uL + uR)); end
                return SVector{1, T}(d / (t * n_sum))
            end
        end
    end
end