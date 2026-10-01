struct LinearAdvection{D, M, T, R} <: HyperbolicPDE{D, M, T, R}
    # Velocity is now explicitly a D-tuple of M x M Jacobians/Speed Matrices
    vel::SVector{D, SMatrix{M, M, T}}
    rep::R
end

# Smart constructor mapping scalars/vectors to proper M x M SMatrices
function LinearAdvection(velocities::Tuple; rep::R = Conservative()) where {R <: EquationRepresentation}

    vel_svec = SVector(velocities...)
    D = length(vel_svec)
    M = size(vel_svec[1], 1)
    T = eltype(vel_svec[1])
    
    return LinearAdvection{D, M, T, R}(vel_svec, rep)
end

@inline prim2cons(::LinearAdvection, U::State) = U
@inline cons2prim(::LinearAdvection, W::State) = W

@inline function flux(eq::LinearAdvection{D, M, T}, U::State{M, T}) where {D, M, T}
    return Flux{D, M, T}(ntuple(d -> eq.vel[d] * U, Val(D)))
end

@inline function max_eigenvalue(eq::LinearAdvection{D, M, T}, U::State{M, T}, d::Int) where {D, M, T}
    # The max wave speed is the maximum eigenvalue of the M x M advection matrix
    # (For M=1, this trivially reduces to the scalar speed)
    return maximum(abs.(eigvals(eq.vel[d])))
end

# API implementation returning the strict M x M matrix
@inline velocity(eq::LinearAdvection, U::State, d::Int) = eq.vel[d]

# Fulfill the core API for kinetic relaxation speeds
@inline function kinetic_wave_speed(eq::LinearAdvection{D, NK, T, R}, d::Int, k::Int) where {D, NK, T, R}
    return eq.vel[d][k,k]
end

function build_equation(::Val{:linear}, pde_conf::Dict, context::Dict)
    rep = parse_representation(pde_conf)
    return LinearAdvection(pde_conf[:velocities]; rep=rep) 
end

# ---------------------------------------------------------
# Analytic Closures
# ---------------------------------------------------------

function analytic_closure(eq::LinearAdvection{D, M, T}, ic::InitialCondition, geom::GeometricDomain) where {D, M, T}
    # 1. Safety Check: Ensure multi-dimensional system matrices commute
    if D > 1
        for i in 1:D
            for j in (i+1):D
                if !isapprox(eq.vel[i] * eq.vel[j], eq.vel[j] * eq.vel[i]; atol=1e-12)
                    error("Velocity matrices for dimensions $i and $j do not commute. Analytic ray-tracing closure is not mathematically valid for this system.")
                end
            end
        end
    end

    mins = geom.mins
    maxs = geom.maxs
    is_per = geom.is_periodic
    
    # 2. System Eigendecomposition (Safe to sum since they commute)
    vel_sum = sum(eq.vel)
    F = eigen(Matrix(vel_sum)) 
    
    # Force the eigen decomposition back into the requested precision T
    R = SMatrix{M, M, T}(real.(F.vectors))
    L = inv(R)
    
    wave_speeds = SVector{M, SVector{D, T}}(ntuple(Val(M)) do m
        SVector{D, T}(ntuple(Val(D)) do d
            (L * eq.vel[d] * R)[m, m]
        end)
    end)
    
    return function exact_linear_system(st::SVector)
        # Ensure spacetime vector inputs are properly cast to T
        t = T(st[end])
        pos = SVector{D, T}(ntuple(d -> T(st[d]), Val(D))) 
        
        u_final = zeros(MVector{M, T})
        
        # 3. Characteristic Tracing
        for m in 1:M
            vel_m = wave_speeds[m]
            
            pos0 = pos - vel_m * t
            
            pos0 = SVector{D, T}(ntuple(Val(D)) do d
                if is_per[d]
                    mins[d] + mod(pos0[d] - mins[d], maxs[d] - mins[d])
                else
                    pos0[d]
                end
            end)
            
            u0 = ic(pos0)
            
            # Project into characteristic variables
            w_m = zero(T)
            for k in 1:M
                w_m += L[m, k] * u0[k]
            end
            
            # Reconstruct into physical variables
            for k in 1:M
                u_final[k] += R[k, m] * w_m
            end
        end
        
        return State{M, T}(Tuple(u_final))
    end
end