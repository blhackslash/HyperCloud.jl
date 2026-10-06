export EquationRepresentation, NCRepresentation, Conservative, analytical_solution

"""
    EquationRepresentation

Abstract base type for categorizing the mathematical formulation of a PDE system.

# Subtypes
- `Conservative`: Concrete type indicating the system evaluates standard conservative fluxes and eigenvalues.
- `NCRepresentation`: Abstract type for non-conservative evaluations (e.g., primitive or Lagrangian representations).
"""
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

include("burgers.jl")
include("linear_advection.jl")
include("euler.jl")