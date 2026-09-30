module HyperCloud

abstract type InitialCondition end
abstract type SmoothInitialCondition <: InitialCondition end

function setInitialConditions!(pg::ParticleGrid{D, M}, eq::HyperbolicPDE, IC::InitialCondition) where {D, M}
    positions = pg.core.positions
    
    @inbounds for i in 1:pg.meta.N
        # Strict assignment: IC(pos) must return an State{M}
        pg.rhos[i] = IC(positions[i])
    end
    return nothing
end

"""
    build_initial_condition(params::ParamDict, D::Int, ::Type{T}) where {T}

Reads the parameters dict and dispatches to the correct IC builder using `Val`.
"""
function build_initial_condition(params::ParamDict, D::Int, ::Type{T}) where {T}
    init_name = lowercase(string(params[:init_func]))
    p = get(params, :init_params, nothing)
    
    # Dispatch on the symbol of the name for modular extensibility
    return build_ic(Val(Symbol(init_name)), p, D, T)
end

# Fallback error for missing implementations
build_ic(name::Val, p, D::Int, ::Type{T}) where {T} = error("Unknown initFunc name: $(typeof(name))")

"""
    generate_analytical_solution(shared_params::ParamDict)

Reads the shared simulation parameters, constructs the corresponding InitialCondition
struct, and returns a fast, standalone closure `exact_u(st)` that evaluates the 
exact analytical solution using a unified spacetime tensor.
"""
"""
    analytical_solution(shared_params::ParamDict)

Reads the shared simulation parameters, constructs the corresponding InitialCondition
and continuous GeometricDomain, and returns a fast, standalone closure `exact_u(st)` 
that evaluates the exact analytical solution using a unified spacetime tensor.
"""
function analytical_solution(shared_params::ParamDict)
    # 1. Instantiate the PDE struct
    # We use Float64 as the standard type for analytical solutions
    eq, D, NM, vel_var = build_equation(shared_params, Float64) 
    
    if !haskey(shared_params, :init_func)
        error("Analytical Factory: Missing 'init_func' in shared_params")
    end
    
    # 2. Build the pure mathematical geometry (handles boundaries and periodicity)
    geom = build_geometric_domain(shared_params, D, Float64)
    
    # 3. Instantiate the InitialCondition struct using the modular builder
    ic = build_initial_condition(shared_params, D, Float64)
    
    # 4. Dispatch to the correct pure mathematical closure 
    # (Notice how it no longer depends on the messy parameters dict)
    return analytic_closure(eq, ic, geom)
end

# ==============================================================================
# PDE-SPECIFIC CLOSURE GENERATORS
# ==============================================================================

# Generic fallback
analytic_closure(eq::HyperbolicPDE, ic::InitialCondition, params::ParamDict) = @warn "No Analytical Solution implemented for $(typeof(eq)) with $(typeof(ic))!"


end