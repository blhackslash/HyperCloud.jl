export build_scheme

# --- Main Scheme Builder ---
"""
    build_scheme(conf::Dict, context::Dict)

Constructs the primary spatial numerical scheme.
Dispatches based on the `:name` key (e.g., `:MUSCL`, `:Upwind`, `:Central`, `:WENO`). Retrieves pre-built dependencies like the selected flux and limiter from the execution context and extracts the scheme's formal `:order` and `:MLS_order` (divergence order).
"""
function build_scheme(conf::Dict, context::Dict)
    name = conf[:name]::Symbol
    return build_scheme(Val(name), conf, context)
end

"""
    build_scheme(name::Val, conf::Dict, context::Dict)

Generic fallback for custom spatial numerical schemes.
Users can extend this by defining `build_scheme(::Val{:my_scheme}, ...)`.
Retrieves pre-built dependencies like the selected flux and limiter from the execution context and extracts the scheme's formal `:order` and `:MLS_order`.
"""
build_scheme(name::Val, conf::Dict, context::Dict) = error("Unknown Scheme: $(typeof(name))")

function build_scheme(::Val{:MUSCL}, conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    
    flux = context[:Flux]
    limiter = context[:Limiter]
    
    order = conf[:order]::Int
    div_order = get(conf,:MLS_order,0)::Int 
    
    return MUSCL(T, D, M, order, limiter, flux; div_order=div_order)
end

function build_scheme(::Val{:Upwind}, conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    
    flux = context[:Flux]
    algType = conf[:upwind_alg_nd]::Symbol 

    order = conf[:order]::Int
    div_order = get(conf,:MLS_order,0)::Int 
    
    return UpwindDivergence(T, D, M, order, algType, flux; div_order=div_order)
end

function build_scheme(::Val{:Central}, conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    
    order = conf[:order]::Int
    div_order = get(conf,:MLS_order,0)::Int 
    
    return CentralDivergence(T, D, M, order; div_order=div_order)
end

function build_scheme(::Val{:WENO}, conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    
    order = conf[:order]::Int
    div_order = get(conf,:MLS_order,0)::Int 
    
    return WENO(T, D, M, order; div_order=div_order)
end

