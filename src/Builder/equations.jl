export build_equation

# --- Main Equation Builder ---
"""
    build_equation(pde_conf::Dict, context::Dict)

Constructs the physical equation system based on the provided configuration dictionary.
Strictly requires a `:name` Symbol (e.g., `:linear`, `:burgers`, `:euler`) in the configuration to dynamically dispatch to the correct specific constructor using `Val`.
"""
function build_equation(pde_conf::Dict, context::Dict)
    if !haskey(pde_conf, :name)
        error("PDE configuration must include a strictly typed :name Symbol (e.g., :linear, :burgers).")
    end
    
    eq_name = pde_conf[:name]::Symbol
    return build_equation(Val(eq_name), pde_conf, context)
end


"""
    build_equation(eq_name::Val, pde_conf::Dict, context::Dict)

Constructs the physical equation system based on the provided configuration dictionary.
Strictly requires a `:name` Symbol (e.g., `:linear`, `:burgers`, `:euler`) in the configuration to dynamically dispatch to the correct specific constructor using `Val`.
"""
build_equation(eq_name::Val, pde_conf::Dict, context::Dict) = error("PDE '$(typeof(eq_name))' is not implemented.")

function build_equation(::Val{:burgers}, pde_conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = pde_conf[:D]::Int
    
    return BurgersEquation{D, T}()
end

function build_equation(::Val{:euler}, pde_conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = pde_conf[:D]::Int
    rep = parse_representation(pde_conf)
    
    gamma = T(pde_conf[:gamma])
    
    return EulerEquation(Val(D), T, gamma, rep)
end

function build_equation(::Val{:linear}, pde_conf::Dict, context::Dict)
    T = context[:T]::DataType
    
    # We pass T securely, and the smart constructor utilizes param2vel internally
    return LinearAdvection(pde_conf[:velocities], T) 
end