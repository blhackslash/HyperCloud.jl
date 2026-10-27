export BGKRelaxation, ExplicitBGKRelaxation

# =========================================================================
# BGK COLLISION TERM
# =========================================================================
"""
    BGKRelaxation(eq::BGKEquation) <: AbstractImplicitSourceTerm
    ExplicitBGKRelaxation(eq::BGKEquation) <: AbstractExplicitSourceTerm

The BGK collision term `S(f) = (M[f] - f) / τ` with the conservative discrete Maxwellian `M[f]` of
[`equilibrium`](@ref) and `τ` from the equation's [`RelaxationTimeModel`](@ref), treated implicitly
(IMEX tableaus) or explicitly (explicit RK tableaus, only stable for `Δt ≲ τ`).

The implicit stage equation `f = f* + Δt a_ii S(f)` is solved in closed form: the collision conserves
the moments, so `M[f] = M[f*]` and `τ(f) = τ(f*)`, and `f = (f* + λ M[f*]) / (1 + λ)` with
`λ = Δt a_ii / τ`. This is exact for any stiffness (`τ → 0` projects onto the equilibrium), which makes
the IMEX scheme asymptotic preserving towards the Euler equations. A non-physical particle state
(`ρ ≤ 0` or `θ ≤ 0`) is left untouched (zero source).
"""
struct BGKRelaxation{E <: BGKEquation} <: AbstractImplicitSourceTerm
    eq::E
end

"The explicit form of the BGK collision term, see [`BGKRelaxation`](@ref)."
struct ExplicitBGKRelaxation{E <: BGKEquation} <: AbstractExplicitSourceTerm
    eq::E
end

const _BGKCollision = Union{BGKRelaxation, ExplicitBGKRelaxation}

@inline function evaluate_source(st::_BGKCollision, f::State{M, T}, p_idx::Int, pg, t::Real) where {M, T}
    U = moments(st.eq, f)
    f_eq, ok = equilibrium(st.eq, U)
    ok || return zero(f)
    return (f_eq - f) / relaxation_time(st.eq, U)
end

@inline function implicit_solve(st::BGKRelaxation, f::State{M, T}, dt_coeff::Real, p_idx::Int, pg, t::Real) where {M, T}
    U = moments(st.eq, f)
    f_eq, ok = equilibrium(st.eq, U)
    ok || return f
    λ = T(dt_coeff) / relaxation_time(st.eq, U)
    return (f + λ * f_eq) / (one(T) + λ)
end

"""
    bgk_sources(eq::BGKEquation, tableau::IMEXButcherTableau) -> (explicit, implicit)

The source tuples of a BGK run: the collision term implicit for an IMEX tableau, explicit for a purely
explicit one, and none in the collisionless limit (`ConstantRelaxationTime(Inf)`).
"""
function bgk_sources(eq::BGKEquation, tableau::IMEXButcherTableau)
    _is_collisionless(eq.relaxation) && return (), ()
    return all(iszero, tableau.a) ? ((ExplicitBGKRelaxation(eq),), ()) : ((), (BGKRelaxation(eq),))
end

# =========================================================================
# INITIAL CONDITIONS
# =========================================================================

"""
    set_initial_conditions!(pg::ParticleGrid, eq::BGKEquation, IC::InitialCondition)

For the BGK equation `IC` gives the macroscopic state `(ρ, ρu..., E)` (the conservative Euler variables
of [`euler_limit`](@ref)); every particle starts in the conservative discrete Maxwellian of that state.
"""
function set_initial_conditions!(pg::ParticleGrid, eq::BGKEquation{D, M, T}, IC::InitialCondition) where {D, M, T}
    set_states!(pg) do pos
        U = IC(pos)
        length(U) == D + 2 || error("BGK initial conditions give the macroscopic state (ρ, ρu..., E) of length $(D + 2), got $(length(U))")
        f, ok = equilibrium(eq, SVector{D + 2, T}(U))
        ok || error("initial state $U at $pos is not admissible for the BGK equation (ρ > 0, θ > 0)")
        f
    end
    return nothing
end
