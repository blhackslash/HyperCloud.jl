export InitialCondition, set_initial_conditions!, build_ic

"""
    InitialCondition

Abstract base type for all initial condition functors in the framework.
Functors derived from this type evaluate the starting state of the physical system at a given spatial coordinate.
"""
abstract type InitialCondition end

"""
    set_initial_conditions!(pg::ParticleGrid, eq::HyperbolicPDE, IC::InitialCondition)
    set_initial_conditions!(pg::ParticleGrid, st::RelaxationSourceTerm, IC::InitialCondition, eq_macro::HyperbolicPDE)

Evaluates the initial condition functor and applies it to the particle grid.

- **Direct formulation:** evaluates `IC` at every particle position and assigns the physical state to `pg.rhos`.
- **Kinetic formulation:** evaluates the macroscopic state and physical fluxes at every particle and initializes
  the kinetic variables directly in the Maxwellian equilibrium, component by component.
"""
function set_initial_conditions!(pg::ParticleGrid{D, M}, eq::HyperbolicPDE, IC::InitialCondition) where {D, M}
    # IC(pos) must return a State{M}
    set_states!(IC, pg)
    return nothing
end

"""
    build_ic(ic_conf::Dict, context::Dict)

Reads the `:IC` namespace and dispatches to the matching builder via `Val(ic_conf[:name])`.
Extend it by defining `build_ic(::Val{:my_ic}, conf::Dict, ctx::Dict)`.
"""
function build_ic(ic_conf::Dict, context::Dict)
    if !haskey(ic_conf, :name)
        error("Initial Condition configuration must include a strictly typed :name Symbol (e.g., :gauss, :box).")
    end

    ic_name = ic_conf[:name]::Symbol
    return build_ic(Val(ic_name), ic_conf, context)
end

# Fallback error for missing implementations
build_ic(name::Val, conf::Dict, ctx::Dict) = error("Unknown initial condition: $(typeof(name))")

include("box.jl")
include("gauss.jl")
include("quadrant_riemann.jl")
include("riemann.jl")
include("sine.jl")
include("smoothed_riemann.jl")
