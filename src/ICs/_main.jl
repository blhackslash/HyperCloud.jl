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

- **Direct Formulation:** Iterates over the grid and evaluates the `IC` functor for each spatial position, assigning the exact physical state to `pg.rhos`[cite: 35].
- **Kinetic Formulation:** Accepts a macroscopic `InitialCondition` and a `RelaxationSourceTerm`. It evaluates the macroscopic state and physical fluxes at each coordinate, and directly initializes the particles into the kinetic Maxwellian equilibrium state component-by-component[cite: 13].
"""
function set_initial_conditions!(pg::ParticleGrid{D, M}, eq::HyperbolicPDE, IC::InitialCondition) where {D, M}
    positions = pg.core.positions
    
    @inbounds for i in 1:pg.meta.N
        # Strict assignment: IC(pos) must return an State{M}
        pg.rhos[i] = IC(positions[i])
    end
    return nothing
end

"""
    build_ic(params::ParamDict, D::Int, ::Type{T}) where {T}

Reads the parameters dict and dispatches to the correct IC builder using `Val`.
"""
function build_ic(ic_conf::Dict, context::Dict)
    if !haskey(ic_conf, :name)
        error("Initial Condition configuration must include a strictly typed :name Symbol (e.g., :gauss, :box).")
    end
    
    ic_name = ic_conf[:name]::Symbol
    return build_ic(Val(ic_name), ic_conf, context)
end

# Fallback error for missing implementations
build_ic(name::Val, p, D::Int, ::Type{T}) where {T} = error("Unknown initFunc name: $(typeof(name))")

include("box.jl")
include("gauss.jl")
include("quadrant_riemann.jl")
include("riemann.jl")
include("sine.jl")
include("smoothed_riemann.jl")