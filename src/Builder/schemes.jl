
# --- Main Scheme Builder ---
function build_scheme(scheme_conf::Dict, context::Dict)
    name = scheme_conf[:name]::Symbol
    return build_scheme(Val(name), scheme_conf, context)
end

build_scheme(name::Val, conf::Dict, context::Dict) = error("Unknown Scheme: $(typeof(name))")

function build_scheme(::Val{:MUSCL}, conf::Dict, context::Dict)
    # 1. Pull base types and dimension dependencies from the shared context
    T = context[:Type]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    
    # 2. Pull pre-built objects from the context
    flux = context[:Flux]
    limiter = context[:Limiter]
    mood = context[:MOOD]
    
    # 3. Pull explicit scheme configurations
    order = conf[:order]::Int
    div_order = conf[:MLS_order]::Int
    
    return MUSCL(T, D, M, order; div_order = div_order, flux = flux, limiter = limiter, mood = mood)
end

function build_scheme(::Val{:Upwind}, conf::Dict, context::Dict)
    T = context[:Type]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    
    flux = context[:Flux]
    order = conf[:order]::Int
    algType = conf[:upwind_alg_nd]::String 
    
    return UpwindDivergence(T, D, M, order; flux = flux, algType = algType)
end

function build_scheme(::Val{:Central}, conf::Dict, context::Dict)
    T = context[:Type]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    order = conf[:order]::Int
    
    return CentralDivergence(T, D, M, order)
end

function build_scheme(::Val{:WENO}, conf::Dict, context::Dict)
    T = context[:Type]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    order = conf[:order]::Int
    
    return WENO(T, D, M, order)
end


# --- Flux Builder ---
function build_flux(flux_conf::Dict, context::Dict)
    name = flux_conf[:name]::Symbol
    return build_flux(Val(name), flux_conf, context)
end

build_flux(name::Val, conf::Dict, context::Dict) = error("Unknown Flux: $(typeof(name))")
build_flux(::Val{:Rusanov}, conf::Dict, context::Dict) = RusanovFlux()
build_flux(::Val{:Upwind}, conf::Dict, context::Dict)  = UpwindFlux()


