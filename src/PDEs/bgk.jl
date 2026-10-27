export BGKEquation, RelaxationTimeModel, ConstantRelaxationTime, ViscosityRelaxationTime
export velocity_grid, num_velocities, velocity_node, moments, equilibrium, relaxation_time, euler_limit

# =========================================================================
# RELAXATION TIME MODELS
# =========================================================================

"""
    RelaxationTimeModel

Relaxation time `τ(ρ, θ)` of the BGK collision operator `(M[f] - f) / τ`, evaluated from the
macroscopic density `ρ` and temperature `θ = p / ρ` of a particle. Subtypes implement
`relaxation_time(model, ρ, θ)`; it runs inside kernels.
"""
abstract type RelaxationTimeModel end

"""
    ConstantRelaxationTime(tau)

Constant relaxation time `τ` (the Knudsen number of the non-dimensional BGK equation). `τ = Inf` is the
collisionless (free transport) limit: the runners drop the collision term and `analytic_closure` returns
the exact solution of the discrete-velocity model.
"""
struct ConstantRelaxationTime{T} <: RelaxationTimeModel
    tau::T
end

"""
    ViscosityRelaxationTime(mu_ref, omega)

`τ = μ(θ) / p` with the power-law viscosity `μ = mu_ref θ^ω` (the Chapman–Enskog viscosity of the BGK
model is `τ p`), i.e. `τ = mu_ref θ^(ω - 1) / ρ`. `ω = 1/2` are hard spheres, `ω = 1` Maxwell molecules.
"""
struct ViscosityRelaxationTime{T} <: RelaxationTimeModel
    mu_ref::T
    omega::T
end

"""
    relaxation_time(model::RelaxationTimeModel, ρ, θ)
    relaxation_time(eq::BGKEquation, U::State)

Relaxation time of the BGK collision term for density `ρ` and temperature `θ`, or for the macroscopic
state `U = (ρ, ρu..., E)` of the equation.
"""
@inline relaxation_time(m::ConstantRelaxationTime, ρ, θ) = m.tau
@inline relaxation_time(m::ViscosityRelaxationTime{T}, ρ, θ) where {T} = m.mu_ref * θ^(m.omega - one(T)) / ρ

_is_collisionless(m::ConstantRelaxationTime) = isinf(m.tau)
_is_collisionless(::RelaxationTimeModel) = false

_convert_relaxation(::Type{T}, m::ConstantRelaxationTime) where {T} = ConstantRelaxationTime{T}(T(m.tau))
_convert_relaxation(::Type{T}, m::ViscosityRelaxationTime) where {T} = ViscosityRelaxationTime{T}(T(m.mu_ref), T(m.omega))

# =========================================================================
# VELOCITY GRIDS
# =========================================================================

"""
    velocity_grid(kind::Symbol, Nv::Int, vmax::Real, ::Type{T} = Float64) -> (nodes, weights)

One-dimensional velocity nodes and quadrature weights (`SVector{Nv, T}` each) of a discrete-velocity
model; the BGK equation uses their tensor product in every spatial direction.

- `:uniform`: midpoint rule on `[-vmax, vmax]`, `v_j = -vmax + (j - 1/2) Δv`, weights `Δv = 2 vmax / Nv`.
  Spectrally accurate for Maxwellians that decay inside the interval; the robust choice for strongly
  non-equilibrium flows.
- `:gauss_hermite`: Gauss–Hermite nodes scaled so that the outermost node is `±vmax`, with the weights
  `w_j e^{x_j^2}` that integrate `f` itself. Exact for polynomials times the Gaussian of that scale, so
  near-equilibrium flows need far fewer nodes (the conservative equilibrium keeps the moments exact
  either way). Weights come from the Christoffel function of the Hermite functions (stable for large `Nv`).
"""
function velocity_grid(kind::Symbol, Nv::Int, vmax::Real, ::Type{T} = Float64) where {T}
    Nv >= 2 || throw(ArgumentError("a velocity grid needs at least 2 nodes, got Nv = $Nv"))
    vmax > 0 || throw(ArgumentError("vmax must be positive, got $vmax"))
    V = Float64(vmax)
    if kind === :uniform
        dv = 2V / Nv
        nodes = [-V + (j - 0.5) * dv for j in 1:Nv]
        weights = fill(dv, Nv)
    elseif kind === :gauss_hermite
        x = eigvals(SymTridiagonal(zeros(Nv), [sqrt(j / 2) for j in 1:(Nv - 1)]))
        # w_j e^{x_j^2} = 1 / Σ_k φ_k(x_j)^2 with the orthonormal Hermite functions φ_k (no over- or underflow)
        W = map(x) do xj
            φ_prev, φ = 0.0, π^(-0.25) * exp(-xj^2 / 2)
            s = φ^2
            for k in 1:(Nv - 1)
                φ_prev, φ = φ, sqrt(2 / k) * xj * φ - sqrt((k - 1) / k) * φ_prev
                s += φ^2
            end
            1 / s
        end
        scale = V / maximum(abs, x)
        nodes = scale .* x
        weights = scale .* W
    else
        throw(ArgumentError("unknown velocity grid :$kind (use :uniform or :gauss_hermite)"))
    end
    return SVector{Nv, T}(T.(nodes)), SVector{Nv, T}(T.(weights))
end

# =========================================================================
# THE EQUATION
# =========================================================================

"""
    BGKEquation{D, M, T, NV1, NV, R} <: HyperbolicPDE{D, M, T}
    BGKEquation(::Val{D}, nodes, weights, gamma, relaxation; center = 0)

Discrete-velocity model of the BGK equation `∂_t f + v·∇f = (M[f] - f) / τ` in `D` space dimensions. The
velocity space is the tensor product of the 1-D `nodes` (with quadrature `weights`, see
[`velocity_grid`](@ref)) in every direction, shifted by `center`: `NV = NV1^D` discrete velocities.

The state of a particle is the vector of distribution values at the nodes, so the number of equations is
the number of velocities. Degrees of freedom beyond the `D` resolved velocity directions (the remaining
velocity directions of a 3-D gas and internal ones, `K = 2 / (γ - 1) - D` in total) are carried by a
second distribution `h` (a Chu reduction, `h = ∫ |ξ|²/2 f dξ` over the unresolved directions): for
`K = 0` (`γ = (D + 2) / D`) the state is `f` (`M = NV`), otherwise `(f, h)` (`M = 2 NV`). Any `γ` thus
needs only a `D`-dimensional velocity grid.

The macroscopic state is `U = (ρ, ρu..., E) = Σ_k w_k (1, v_k, |v_k|²/2) f_k + (0, 0, h_k)` (see
[`moments`](@ref)), the conservative variables of the Euler equations with the same `γ`
([`euler_limit`](@ref)), which the BGK solution approaches for `τ → 0`. The collision term relaxes
towards the conservative discrete Maxwellian of [`equilibrium`](@ref) (source terms `BGKRelaxation`,
`ExplicitBGKRelaxation`) with the relaxation time of `relaxation` (a [`RelaxationTimeModel`](@ref)).

The transport is linear and diagonal (velocity `k` moves with `v_k`): the flux is `v_{k,d} f_k`, the
Rusanov and Upwind fluxes reduce to exact upwinding per velocity (dissipation `|v_{k,d}|` instead of the
largest speed for every component), and `max_eigenvalue` is the largest node speed. The equation object
is a small bits type (nodes, weights, center, `K`, `γ`, the relaxation model); node indices are resolved
at compile time, so there are no per-velocity arrays in the kernels.
"""
struct BGKEquation{D, M, T, NV1, NV, R <: RelaxationTimeModel} <: HyperbolicPDE{D, M, T}
    nodes::SVector{NV1, T}
    weights::SVector{NV1, T}
    center::SVector{D, T}
    max_speed::SVector{D, T}
    K::T
    gamma::T
    relaxation::R
end

function BGKEquation(::Val{D}, nodes::SVector{NV1, T}, weights::SVector{NV1, T}, gamma::Real,
                     relaxation::RelaxationTimeModel; center = ntuple(_ -> 0, D)) where {D, NV1, T}
    K = 2 / (Float64(gamma) - 1) - D
    abs(K) <= 1e-8 && (K = 0.0)
    K >= 0 || throw(ArgumentError("γ = $gamma is too large for D = $D (BGK needs γ ≤ (D + 2) / D = $((D + 2) / D))"))
    NV = NV1^D
    M  = iszero(K) ? NV : 2NV
    c  = SVector{D, T}(ntuple(d -> T(center[d]), Val(D)))
    smax = SVector{D, T}(ntuple(d -> maximum(abs, c[d] .+ nodes), Val(D)))
    rel = _convert_relaxation(T, relaxation)
    return BGKEquation{D, M, T, NV1, NV, typeof(rel)}(nodes, weights, c, smax, T(K), T((D + K + 2) / (D + K)), rel)
end

"""
    num_velocities(eq::BGKEquation) -> Int

Number of discrete velocities `NV = NV1^D` (the number of equations is `NV`, or `2 NV` with internal
degrees of freedom).
"""
num_velocities(::BGKEquation{D, M, T, NV1, NV}) where {D, M, T, NV1, NV} = NV

@inline _has_h(::BGKEquation{D, M, T, NV1, NV}) where {D, M, T, NV1, NV} = M != NV

# 1-D node index along `d` of state component `k` (the `h` components share the nodes of `f`), the first
# axis running fastest.
_node_index(NV1::Int, NV::Int, k::Int, d::Int) = ((k - 1) % NV) ÷ NV1^(d - 1) % NV1 + 1

# The node indices of all components as a compile-time constant (a literal tuple): the index arithmetic
# never runs in the kernels, a runtime `k` is a table load and a constant one folds away.
@generated function _node_table(::BGKEquation{D, M, T, NV1, NV}) where {D, M, T, NV1, NV}
    table = ntuple(k -> ntuple(d -> _node_index(NV1, NV, k, d), D), M)
    return :($table)
end
@inline _node_index(eq::BGKEquation, k::Int, d::Int) = @inbounds _node_table(eq)[k][d]

@inline _velocity_component(eq::BGKEquation, k::Int, d::Int) = @inbounds eq.center[d] + eq.nodes[_node_index(eq, k, d)]

"""
    velocity_node(eq::BGKEquation, k::Int) -> (v::Space{D, T}, w::T)

Velocity and quadrature weight of state component `k` (components `NV + 1 … 2 NV` of `h` share the
nodes of `f`).
"""
@inline function velocity_node(eq::BGKEquation{D, M, T}, k::Int) where {D, M, T}
    v = Space{D, T}(ntuple(d -> _velocity_component(eq, k, d), Val(D)))
    w = prod(ntuple(d -> (@inbounds eq.weights[_node_index(eq, k, d)]), Val(D)))
    return v, w
end

# Velocities of all state components along axis `d`, unrolled with literal node indices.
@generated function _component_speeds(eq::BGKEquation{D, M, T, NV1, NV}, ::Val{d}) where {D, M, T, NV1, NV, d}
    terms = [:(eq.center[$d] + eq.nodes[$(_node_index(NV1, NV, k, d))]) for k in 1:M]
    return :(@inbounds State{$M, $T}(($(terms...),)))
end

# `Flux{D, M, T}((g(speeds_1), …, g(speeds_D)))` with the speeds of every axis (`g(speeds, d)`), unrolled over D.
@generated function _per_axis(g::G, eq::BGKEquation{D, M, T}) where {G, D, M, T}
    terms = [:(g(_component_speeds(eq, Val($d)), $d)) for d in 1:D]
    return :(Flux{$D, $M, $T}(($(terms...),)))
end

# =========================================================================
# MOMENTS AND EQUILIBRIUM
# =========================================================================

"""
    moments(eq::BGKEquation, f::State) -> State{D + 2}

Macroscopic conservative state `(ρ, ρu..., E)` of the distribution `f` by quadrature over the velocity
grid (`E` includes the energy `Σ w_k h_k` of the unresolved degrees of freedom).
"""
@inline function moments(eq::BGKEquation{D, M, T, NV1, NV}, f::State{M, T}) where {D, M, T, NV1, NV}
    m = zero(SVector{D + 2, T})
    @inbounds for k in 1:NV
        v, w = velocity_node(eq, k)
        fk = f[k]
        e  = T(0.5) * sum(abs2, v) * fk
        _has_h(eq) && (e += f[k + NV])
        m += w * SVector{D + 2, T}(fk, (fk * v)..., e)
    end
    return m
end

# Density, velocity and temperature θ = p / ρ of a macroscopic state.
@inline function _primitives(eq::BGKEquation{D, M, T}, U::SVector{P, T}) where {D, M, T, P}
    ρ = U[1]
    u = Space{D, T}(ntuple(d -> U[1 + d], Val(D))) / ρ
    θ = T(2) * (U[P] - T(0.5) * ρ * sum(abs2, u)) / ((T(D) + eq.K) * ρ)
    return ρ, u, θ
end

@inline function relaxation_time(eq::BGKEquation{D, M, T}, U::SVector{P, T}) where {D, M, T, P}
    ρ, _, θ = _primitives(eq, U)
    return relaxation_time(eq.relaxation, ρ, θ)
end

# The conservative discrete Maxwellian is f_k = exp(α·ψ_k) (2π θ̂)^(K/2), h_k = (K/2) θ̂ f_k with
# ψ_k = (1, v_k, |v_k|²/2) and θ̂ = -1/α_E. Its moments are the gradient of the convex dual potential
# Φ(α) = Σ_k w_k exp(α·ψ_k) (2π θ̂)^(K/2), so the Newton matrix (the Hessian of Φ) is symmetric positive
# definite. Returns the moments m(α), the Hessian J(α) and the moment scale.
@inline function _equilibrium_system(eq::BGKEquation{D, M, T, NV1, NV}, α::SVector{P, T}) where {D, M, T, NV1, NV, P}
    θh = -one(T) / α[P]
    c  = T(0.5) * eq.K
    lognorm = c * log(T(2π) * θh)
    S0 = zero(T)
    S1 = zero(SVector{P, T})
    S2 = zero(SMatrix{P, P, T})
    @inbounds for k in 1:NV
        v, w = velocity_node(eq, k)
        ψ  = SVector{P, T}(one(T), v..., T(0.5) * sum(abs2, v))
        G  = w * exp(dot(α, ψ) + lognorm)
        S0 += G
        S1 += G * ψ
        S2 += G * (ψ * ψ')
    end
    _has_h(eq) || return S1, S2
    eP = SVector{P, T}(ntuple(i -> i == P ? one(T) : zero(T), Val(P)))
    m  = S1 + (c * θh * S0) * eP
    J  = S2 + (c * θh) * (S1 * eP' + eP * S1') + ((c^2 + c) * θh^2 * S0) * (eP * eP')
    return m, J
end

# Cholesky solve of a small SPD system without allocation (kernel safe); `false` if not positive definite.
@inline function _spd_solve(A::SMatrix{P, P, T}, b::SVector{P, T}) where {P, T}
    L = zero(SMatrix{P, P, T})
    @inbounds for j in 1:P
        s = A[j, j]
        for k in 1:(j - 1)
            s -= L[j, k]^2
        end
        s > zero(T) || return b, false
        ljj = sqrt(s)
        L = setindex(L, ljj, j, j)
        for i in (j + 1):P
            t = A[i, j]
            for k in 1:(j - 1)
                t -= L[i, k] * L[j, k]
            end
            L = setindex(L, t / ljj, i, j)
        end
    end
    y = zero(SVector{P, T})
    @inbounds for i in 1:P
        t = b[i]
        for k in 1:(i - 1)
            t -= L[i, k] * y[k]
        end
        y = setindex(y, t / L[i, i], i)
    end
    x = zero(SVector{P, T})
    @inbounds for i in P:-1:1
        t = y[i]
        for k in (i + 1):P
            t -= L[k, i] * x[k]
        end
        x = setindex(x, t / L[i, i], i)
    end
    return x, true
end

# Moment residual relative to the natural scale of each component (mass ρ, momentum ρ(|u| + √θ), energy E).
@inline function _scaled_residual(r::SVector{P, T}, scale::SVector{P, T}) where {P, T}
    return maximum(abs.(r) ./ scale)
end

const _BGK_NEWTON_ITERATIONS = 12
const _BGK_LINE_SEARCH_STEPS = 6

"""
    equilibrium(eq::BGKEquation, U::State{D + 2}) -> (f_eq::State{M}, ok::Bool)

Conservative discrete Maxwellian of the macroscopic state `U = (ρ, ρu..., E)`: the discrete entropic
equilibrium `f_k = exp(α·ψ_k)` (Mieussens), whose discrete moments reproduce `U` to round-off on the
velocity grid, so the collision term conserves mass, momentum and energy exactly. `α` is found by a
damped Newton iteration on the convex dual problem, started from the continuous Maxwellian (on a grid
that resolves the Maxwellian that start is already accurate to the quadrature error, so one or two
iterations suffice). With internal degrees of freedom `h_k = (K/2) θ̂ f_k`.

`ok` is `false` for a non-physical `U` (non-finite, `ρ ≤ 0` or `θ ≤ 0`); `f_eq` is meaningless then.
If the grid cannot represent `U` (e.g. a mean velocity near `vmax`), the iteration stops at the closest
equilibrium it reaches. Allocation free, runs inside kernels.
"""
function equilibrium(eq::BGKEquation{D, M, T, NV1, NV}, U::SVector{P, T}) where {D, M, T, NV1, NV, P}
    ρ, u, θ = _primitives(eq, U)
    (isfinite(ρ) && isfinite(θ) && all(isfinite, u) && ρ > zero(T) && θ > zero(T)) || return zero(State{M, T}), false

    # Start: the continuous Maxwellian ρ (2πθ)^(-(D+K)/2) exp(-|v - u|²/(2θ)) in the exponential form
    α = SVector{P, T}(log(ρ) - T(0.5) * (T(D) + eq.K) * log(T(2π) * θ) - T(0.5) * sum(abs2, u) / θ,
                      (u / θ)..., -one(T) / θ)
    vscale = sqrt(sum(abs2, u)) + sqrt(θ)
    scale  = SVector{P, T}(ρ, ntuple(_ -> ρ * vscale, Val(D))..., U[P])
    tol    = T(16) * eps(T)

    m, J = _equilibrium_system(eq, α)
    res  = _scaled_residual(m - U, scale)
    for _ in 1:_BGK_NEWTON_ITERATIONS
        res <= tol && break
        Δ, spd = _spd_solve(J, m - U)
        spd || break
        step = one(T)
        accepted = false
        for _ in 1:_BGK_LINE_SEARCH_STEPS
            αn = α - step * Δ
            if all(isfinite, αn) && αn[P] < zero(T)
                mn, Jn = _equilibrium_system(eq, αn)
                resn = _scaled_residual(mn - U, scale)
                if resn < res
                    α, m, J, res = αn, mn, Jn, resn
                    accepted = true
                    break
                end
            end
            step *= T(0.5)
        end
        accepted || break        # round-off floor (or an unrepresentable state): keep the best iterate
    end
    return _equilibrium_state(eq, α), true
end

@inline function _equilibrium_state(eq::BGKEquation{D, M, T, NV1, NV}, α::SVector{P, T}) where {D, M, T, NV1, NV, P}
    θh = -one(T) / α[P]
    c  = T(0.5) * eq.K
    lognorm = c * log(T(2π) * θh)
    return State{M, T}(ntuple(Val(M)) do k
        v, _ = velocity_node(eq, k)
        G = exp(dot(α, SVector{P, T}(one(T), v..., T(0.5) * sum(abs2, v))) + lognorm)
        k > NV ? c * θh * G : G
    end)
end

"""
    euler_limit(eq::BGKEquation) -> EulerEquation

The compressible Euler equations (conservative variables) with the ratio of specific heats of the BGK
model, which its moments approach for `τ → 0`: the reference for small relaxation times.
"""
euler_limit(eq::BGKEquation{D, M, T}) where {D, M, T} = EulerEquation(Val(D), T, eq.gamma, Conservative())

# =========================================================================
# HYPERCLOUDCORE API
# =========================================================================

@inline flux(eq::BGKEquation{D, M, T}, U::State{M, T}) where {D, M, T} = _per_axis((v, d) -> v .* U, eq)

@inline max_eigenvalue(eq::BGKEquation, U::State, d::Int) = @inbounds eq.max_speed[d]

"Admissible BGK state (MOOD PAD): finite, with positive density and temperature (individual `f_k` may dip below zero)."
@inline function is_admissible(eq::BGKEquation{D, M, T}, U::State{M, T}) where {D, M, T}
    all(isfinite, U) || return false
    ρ, _, θ = _primitives(eq, moments(eq, U))
    return ρ > zero(T) && θ > zero(T)
end

@inline fluid_velocity(eq::BGKEquation{D, M, T}, U::State{M, T}) where {D, M, T} = _primitives(eq, moments(eq, U))[2]

# The transport is diagonal, so the Rusanov flux of every velocity is its exact upwind flux: dissipation
# |v_{k,d}| per component, not the largest node speed for all of them (which smears the slow velocities and
# dominates the error near equilibrium). `UpwindFlux` and `AlignedRusanovFlux` fall back to this method.
@inline function (::RusanovFlux)(f_L::Flux{D, M, T}, f_R::Flux{D, M, T}, F_L::Flux{D, M, T}, F_R::Flux{D, M, T},
                                 eq::BGKEquation{D, M, T}) where {D, M, T}
    return _per_axis(eq) do v, d
        T(0.5) * (F_L[d] + F_R[d] - abs.(v) .* (f_R[d] - f_L[d]))
    end
end

# =========================================================================
# ANALYTIC CLOSURE (collisionless limit)
# =========================================================================

"""
    analytic_closure(eq::BGKEquation, ic::InitialCondition, geom::GeometricDomain)

For the collisionless model (`ConstantRelaxationTime(Inf)`) the exact solution of the discrete-velocity
equation: every velocity is transported, `f_k(x, t) = f_k^0(x - v_k t)` with the equilibrium initial data
of `ic` (wrapped on periodic axes), and the closure returns its moments `(ρ, ρu..., E)` like the BGK
runner. Finite relaxation times have no closed-form solution (warns, returns `nothing`); compare with the
Euler limit `analytic_closure(euler_limit(eq), ic, geom)` for small `τ` instead.
"""
function analytic_closure(eq::BGKEquation{D, M, T, NV1, NV}, ic::InitialCondition, geom::GeometricDomain) where {D, M, T, NV1, NV}
    if !_is_collisionless(eq.relaxation)
        @warn "No Analytical Solution implemented for the BGK equation with a finite relaxation time; " *
              "compare with analytic_closure(euler_limit(eq), ic, geom) for small τ."
        return nothing
    end
    mins, maxs, is_per = geom.mins, geom.maxs, geom.is_periodic

    return function exact_free_transport(st::SVector)
        t   = T(st[end])
        pos = Space{D, T}(ntuple(d -> T(st[d]), Val(D)))
        U   = zero(SVector{D + 2, T})
        for k in 1:NV
            v, w = velocity_node(eq, k)
            x0 = pos - v * t
            x0 = Space{D, T}(ntuple(Val(D)) do d
                is_per[d] ? T(mins[d] + mod(x0[d] - mins[d], maxs[d] - mins[d])) : x0[d]
            end)
            f0, ok = equilibrium(eq, SVector{D + 2, T}(ic(x0)))
            ok || error("initial state at $x0 is not admissible for the BGK equation")
            e = T(0.5) * sum(abs2, v) * f0[k]
            _has_h(eq) && (e += f0[k + NV])
            U += w * SVector{D + 2, T}(f0[k], (f0[k] * v)..., e)
        end
        return U
    end
end
