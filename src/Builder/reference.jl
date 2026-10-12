export analytic_closure, analytical_solution

"""
    analytical_solution(shared_params::ParamDict)

Builds the equation, domain and initial condition from the `PDE`, `Grid` and `IC`
namespaces and returns the standalone closure `exact_u(st)` from [`analytic_closure`](@ref).
`st` is the spacetime vector `SVector(x..., t)`.

The reference is computed in `shared_params[:real_type]` (default `Float64`), independently
of the simulation precision `:T`, so Float32 runs are still compared to a Float64 reference.
"""
function analytical_solution(shared_params::ParamDict)
    # 1. Precision of the reference (independent of the simulation precision :T)
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

Returns a closure `st -> State{M}` evaluating the exact solution at the spacetime point
`st = SVector(x..., t)`. PDEs add methods for the initial conditions they can solve exactly.

This generic fallback only emits a warning and returns `nothing` when no method exists for
the combination of equation, initial condition and domain.
"""
analytic_closure(eq::HyperbolicPDE, ic::InitialCondition, geom::GeometricDomain) = @warn "No Analytical Solution implemented for $(typeof(eq)) with $(typeof(ic)) on domain $(typeof(geom))!"
