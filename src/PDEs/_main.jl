export EquationRepresentation, NCRepresentation, Conservative

abstract type EquationRepresentation end
abstract type NCRepresentation <: EquationRepresentation end
struct Conservative <: EquationRepresentation end

# --- Representation Parsing ---
function parse_representation(pde_conf::Dict)
    rep_sym = get(pde_conf, :representation, :conservative)::Symbol
    return _parse_representation(Val(rep_sym), pde_conf)
end

_parse_representation(::Val{:conservative}, pde_conf::Dict) = Conservative()
function _parse_representation(::Val{:primitive}, pde_conf::Dict)
    path_sym = get(pde_conf, :path, :mapped)::Symbol
    return Primitive(_parse_path(Val(path_sym)))
end
function _parse_representation(::Val{:lagrangian}, pde_conf::Dict)
    path_sym = get(pde_conf, :path, :mapped)::Symbol
    return Lagrangian(_parse_path(Val(path_sym)))
end
_parse_representation(rep_val::Val, pde_conf::Dict) = error("Unknown PDE representation: $(typeof(rep_val))")

_parse_path(::Val{:line}) = LinePath()
_parse_path(::Val{:mapped}) = MappedPath()
_parse_path(::Val{:naive}) = NaiveAveragePath()
_parse_path(path_val::Val) = error("Unknown PDE path: $(typeof(path_val))")

# --- Main Equation Builder ---
function build_equation(pde_conf::Dict, context::Dict)
    if !haskey(pde_conf, :name)
        error("PDE configuration must include a strictly typed :name Symbol (e.g., :linear, :burgers).")
    end
    
    eq_name = pde_conf[:name]::Symbol
    return build_equation(Val(eq_name), pde_conf, context)
end

# Generic fallback
build_equation(eq_name::Val, pde_conf::Dict, context::Dict) = error("PDE '$(typeof(eq_name))' is not implemented.")

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
    context[:Type] = T

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
analytic_closure(eq::HyperbolicPDE, ic::InitialCondition, params::ParamDict) = @warn "No Analytical Solution implemented for $(typeof(eq)) with $(typeof(ic))!"

include("burgers.jl")
include("linear_advection.jl")
include("euler.jl")