export analytic_closure, analytical_solution

"""
    analytical_solution(shared_params::ParamDict)

Reads the shared simulation parameters, constructs the corresponding InitialCondition
and continuous GeometricDomain, and returns a fast, standalone closure `exact_u(st)` 
that evaluates the exact analytical solution using a unified spacetime tensor.
"""
function analytical_solution(shared_params::ParamDict)
    # 1. Instantiate the PDE struct
    # We use Float64 as the standard type for analytical solutions
    T = get(shared_params, :real_type, Float64) 
    context = Dict{Symbol, Any}()
    context[:T] = T

    pde_conf    = extract_namespace(shared_params, :PDE)
    domain_conf = extract_namespace(shared_params, :Grid)
    ic_conf     = extract_namespace(shared_params, :IC)

    # 2. Base Equation & Dimensions
    eq = build_equation(pde_conf, context)
    D = get_D(eq)
    
    context[:Equation] = eq
    context[:D] = D
    context[:M] = get_M(eq)
    
    # 2. Build the pure mathematical geometry (handles boundaries and periodicity)
    geom = build_domain(domain_conf, context)
    context[:Domain] = geom
    
    # 3. Instantiate the InitialCondition struct using the modular builder
    IC = build_ic(ic_conf, context)
    
    # 4. Dispatch to the correct pure mathematical closure 
    # (Notice how it no longer depends on the messy parameters dict)
    return analytic_closure(eq, IC, geom)
end

# ==============================================================================
# PDE-SPECIFIC CLOSURE GENERATORS
# ==============================================================================

# Generic fallback
"""
    analytic_closure(eq::HyperbolicPDE, ic::InitialCondition, geom::GeometricDomain)

Generates an exact mathematical closure evaluating the analytical state for a given PDE, initial condition, and continuous geometric domain.
This function serves as a generic fallback that emits a warning if no exact analytical solution has been explicitly implemented for the provided combination of physics, geometry, and initial states.
"""
analytic_closure(eq::HyperbolicPDE, ic::InitialCondition, geom::GeometricDomain) = @warn "No Analytical Solution implemented for $(typeof(eq)) with $(typeof(ic)) on domain $(typeof(geom))!"
