export build_domain
# =========================================================================
# MODULAR DOMAIN BUILDER
# =========================================================================

# --- Boundary Condition Parsing ---
function parse_bc(domain_conf::Dict)
    raw_bc = get(domain_conf, :bc, Dict{Int, AbstractBoundaryCondition}())
    
    if raw_bc isa String
        raw_bc = eval(Meta.parse(raw_bc))
    end
    
    bc_map = Dict{Int, AbstractBoundaryCondition}()
    
    for (tag, bc_obj) in raw_bc
        if bc_obj isa AbstractBoundaryCondition
            bc_map[tag] = bc_obj
        else
            bc_sym = Symbol(bc_obj)
            if bc_sym === :outflow
                bc_map[tag] = OutflowBC()
            elseif bc_sym === :fixed_dirichlet
                bc_map[tag] = FixedDirichlet()
            else
                if isdefined(Main, bc_sym)
                    bc_map[tag] = getfield(Main, bc_sym)()
                elseif isdefined(@__MODULE__, bc_sym)
                    bc_map[tag] = getfield(@__MODULE__, bc_sym)()
                else
                    error("Boundary Condition '$bc_sym' could not be found.")
                end
            end
        end
    end
    
    return bc_map
end

# --- Main Domain Builder ---
"""
    build_domain(domain_conf::Dict, context::Dict)

Instantiates the continuous `GeometricDomain` (e.g., `:rectangular` or `:spherical`) and parses boundary conditions.
Requires a `:domain` key unless a pre-instantiated `GeometricDomain` is passed via the `:instance` key. Dynamically resolves boundary conditions using `parse_bc` and maps periodic boundary flags.
"""
function build_domain(domain_conf::Dict, context::Dict)
    # Direct pass-through if the user supplied an already-instantiated GeometricDomain
    if domain_conf[:domain] isa GeometricDomain
        return domain_conf[:domain]
    end
    
    if !haskey(domain_conf, :domain)
        error("Domain configuration must strictly include a :domain key (e.g., :rectangular, :spherical).")
    end
    
    domain_name = domain_conf[:domain]::Symbol
    return build_domain(Val(domain_name), domain_conf, context)
end

"""
    build_domain(name::Val, domain_conf::Dict, context::Dict)

Generic fallback for custom domain instantiation. 
Users can extend this by defining `build_domain(::Val{:my_domain}, ...)`. 
Instantiates the continuous `GeometricDomain` (e.g., `:rectangular` or `:spherical`) and parses boundary conditions.
"""
build_domain(name::Val, domain_conf::Dict, context::Dict) = error("Unknown domain shape: $(typeof(name))")

# --- Specific Domain Builders ---
function build_domain(::Val{:rectangular}, domain_conf::Dict, context::Dict)
    T = context[:T]::DataType
    
    bc_map = parse_bc(domain_conf)
    
    # Directly pull periodicity (expected to be Bool or NTuple{D, Bool})
    is_per = get(domain_conf, :periodic, false)
    
    req_mins = Tuple(T.(domain_conf[:mins]))
    req_maxs = Tuple(T.(domain_conf[:maxs]))
    
    return get_rectangular_domain(T, req_mins, req_maxs; bc_map = bc_map, is_periodic = is_per)
end

function build_domain(::Val{:spherical}, domain_conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = context[:D]::Int
    
    bc_map = parse_bc(domain_conf)
    
    # Directly pull periodicity
    is_per = get(domain_conf, :periodic, false)
    
    req_mins = Tuple(T.(domain_conf[:mins]))
    req_maxs = Tuple(T.(domain_conf[:maxs]))
    
    center = ntuple(d -> (req_mins[d] + req_maxs[d]) / 2.0, Val(D))
    radius = (req_maxs[1] - req_mins[1]) / 2.0
    
    return get_spherical_domain(T, center, radius; bc_map = bc_map, is_periodic = is_per)
end