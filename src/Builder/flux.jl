export build_flux

# --- Flux Builder ---
"""
    build_flux(flux_conf::Dict, context::Dict)

Instantiates the numerical flux function. 
Supports options like `:Rusanov` or `:Upwind` based on the `:name` key.
"""
function build_flux(flux_conf::Dict, context::Dict)
    name = flux_conf[:name]::Symbol
    return build_flux(Val(name), flux_conf, context)
end

"""
    build_flux(name::Val, conf::Dict, context::Dict)

Generic fallback for custom numerical flux functions. 
Users can extend this by defining `build_flux(::Val{:my_flux}, ...)`.
"""
build_flux(name::Val, conf::Dict, context::Dict) = error("Unknown Flux: $(typeof(name))")
build_flux(::Val{:Rusanov}, conf::Dict, context::Dict) = RusanovFlux()
build_flux(::Val{:Upwind}, conf::Dict, context::Dict)  = UpwindFlux()