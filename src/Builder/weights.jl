export build_weights

# =========================================================================
# MODULAR WEIGHT BUILDER
# =========================================================================
"""
    build_weights(weight_conf::Dict, context::Dict)

Constructs the meshfree Moving Least Squares (MLS) weight function.
Requires a `:name` (e.g., `:exponential` or `:inverse`), an interpolation `:range` scalar, and an `:alpha` parameter. Computes the absolute interpolation range dynamically by multiplying the requested `:range` with the `max_dx` extracted from the execution context.
"""
function build_weights(weight_conf::Dict, context::Dict)
    if !haskey(weight_conf, :name)
        error("Weight configuration must include a strictly typed :name Symbol (e.g., :exponential, :inverse).")
    end
    if !haskey(weight_conf, :range)
        error("Weight configuration must include a strictly typed :range parameter for the interpolation scaling factor.")
    end
    
    weight_name = weight_conf[:name]::Symbol
    
    # Extract dependencies from the context
    T = context[:T]::DataType
    max_dx = context[:max_dx]::T
    
    # Compute the absolute interpolation range directly in the weight builder
    # and save it to the context for standard API retrieval
    context[:interp_range] = T(weight_conf[:range]) * max_dx
    
    return build_weights(Val(weight_name), weight_conf, context)
end

# Generic fallback
"""
    build_weights(name::Val, weight_conf::Dict, context::Dict)

Generic fallback for custom meshfree MLS weight functions.
Users can extend this by defining `build_weights(::Val{:my_weight}, ...)`.
Requires a `:range` scalar and an `:alpha` parameter, computing the absolute interpolation range dynamically by multiplying the requested `:range` with `max_dx` from the execution context.
"""
build_weights(name::Val, weight_conf::Dict, context::Dict) = error("Unknown weight function: $(typeof(name))")

# --- Specific Weight Builders ---
function build_weights(::Val{:exponential}, weight_conf::Dict, context::Dict)
    T = context[:T]::DataType
    interp_range = context[:interp_range]::T
    
    # Explicitly require alpha (no defaults!)
    alpha = T(weight_conf[:alpha])
    return ExponentialWeightFunction(alpha, interp_range)
end

function build_weights(::Val{:inverse}, weight_conf::Dict, context::Dict)
    T = context[:T]::DataType
    interp_range = context[:interp_range]::T
    
    # Explicitly require alpha (no defaults!)
    alpha = T(weight_conf[:alpha])
    return InverseWeightFunction(alpha, interp_range)
end