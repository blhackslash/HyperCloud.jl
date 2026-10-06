export build_limiter

# --- Limiter Builder ---
"""
    build_limiter(limiter_conf::Dict, context::Dict)

Builds the spatial slope limiter based on the `:name` key.
Supports `:none`, `:minmod`, `:superbee`, `:VK` (Venkatakrishnan), and `:BJ` (Barth-Jespersen) limiters. All active limiters strictly require a `:mode` parameter (e.g., `:hard`).
"""
function build_limiter(limiter_conf::Dict, context::Dict)
    name = get(limiter_conf,:name,:none)::Symbol
    return build_limiter(Val(name), limiter_conf, context)
end
"""
    build_limiter(name::Val, conf::Dict, context::Dict)

Generic fallback for custom spatial slope limiters.
Users can extend this by defining `build_limiter(::Val{:my_limiter}, ...)`.
Supports built-in options like `:none`, `:minmod`, `:superbee`, `:VK`, and `:BJ`, which strictly require a `:mode` parameter (e.g., `:hard`).
"""
build_limiter(name::Val, conf::Dict, context::Dict) = error("Unknown Limiter: $(typeof(name))")
build_limiter(::Val{:none}, conf::Dict, context::Dict) = NoLimiter()

function build_limiter(::Val{:minmod}, conf::Dict, context::Dict)
    return MinmodLimiter(conf[:mode]::Symbol)
end
function build_limiter(::Val{:superbee}, conf::Dict, context::Dict)
    return SuperbeeLimiter(conf[:mode]::Symbol)
end
function build_limiter(::Val{:VK}, conf::Dict, context::Dict)
    return VenkatakrishnanLimiter(conf[:mode]::Symbol)
end
function build_limiter(::Val{:BJ}, conf::Dict, context::Dict)
    return BarthJespersenLimiter(conf[:mode]::Symbol)
end
