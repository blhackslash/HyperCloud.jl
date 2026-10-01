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
            if bc_sym === :outflow || bc_sym === :OutflowBC
                bc_map[tag] = OutflowBC()
            elseif bc_sym === :fixed_dirichlet || bc_sym === :FixedDirichlet
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

# --- Periodicity Parsing ---
function parse_periodic(domain_conf::Dict, D::Int)
    is_per_input = get(domain_conf, :periodic, false)
    return isa(is_per_input, Bool) ? SVector{D, Bool}(ntuple(_ -> is_per_input, Val(D))) : SVector{D, Bool}(is_per_input)
end

# --- Main Domain Builder ---
function build_domain(domain_conf::Dict, D::Int, ::Type{T}) where {T}
    # Direct pass-through if the user supplied an already-instantiated GeometricDomain
    if haskey(domain_conf, :instance) && domain_conf[:instance] isa GeometricDomain
        return domain_conf[:instance]
    end
    
    domain_name = get(domain_conf, :name, :rectangular)::Symbol
    return build_domain(Val(domain_name), domain_conf, D, T)
end

# Generic fallback
build_domain(name::Val, domain_conf::Dict, D::Int, ::Type{T}) where {T} = error("Unknown domain shape: $(typeof(name))")

# --- Specific Domain Builders ---
function build_domain(::Val{:rectangular}, domain_conf::Dict, D::Int, ::Type{T}) where {T}
    bc_map = parse_bc(domain_conf)
    is_per = parse_periodic(domain_conf, D)
    
    req_mins = T.(domain_conf[:mins]::Tuple)
    req_maxs = T.(domain_conf[:maxs]::Tuple)
    
    return get_rectangular_domain(T, req_mins, req_maxs; bc_map = bc_map, is_periodic = is_per)
end

function build_domain(::Val{:spherical}, domain_conf::Dict, D::Int, ::Type{T}) where {T}
    bc_map = parse_bc(domain_conf)
    is_per = parse_periodic(domain_conf, D)
    
    req_mins = T.(domain_conf[:mins]::Tuple)
    req_maxs = T.(domain_conf[:maxs]::Tuple)
    
    center = ntuple(d -> (req_mins[d] + req_maxs[d]) / 2.0, Val(D))
    radius = (req_maxs[1] - req_mins[1]) / 2.0
    
    return get_spherical_domain(T, center, radius; bc_map = bc_map, is_periodic = is_per)
end