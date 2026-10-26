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

"""
    build_equation(::Val{:bgk}, pde_conf::Dict, context::Dict) -> BGKEquation

Discrete-velocity BGK model. Keys of the `PDE` namespace:

| Key | Meaning |
|:--|:--|
| `:D` | space dimension (the velocity grid has the same dimension) |
| `:Nv` | velocity nodes per direction (`Nv^D` discrete velocities) |
| `:vmax` | half-width of the velocity grid (outermost node for `:gauss_hermite`) |
| `:velocity_grid` | `:uniform` (default) or `:gauss_hermite`, see [`velocity_grid`](@ref) |
| `:vcenter` | centre of the velocity grid per direction (default `0`) |
| `:gamma` | ratio of specific heats (default `(D + 2) / D`, no internal degrees of freedom) |
| `:tau` | constant relaxation time ([`ConstantRelaxationTime`](@ref)); `Inf` is collisionless |
| `:viscosity`, `:omega` | instead of `:tau`: `τ = μ / p` with `μ = viscosity θ^omega` ([`ViscosityRelaxationTime`](@ref), `omega` defaults to `0.5`) |
"""
function build_equation(::Val{:bgk}, pde_conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = pde_conf[:D]::Int
    haskey(pde_conf, :tau) && haskey(pde_conf, :viscosity) && error("BGK: give either :tau or :viscosity, not both.")
    relaxation = if haskey(pde_conf, :viscosity)
        ViscosityRelaxationTime(T(pde_conf[:viscosity]), T(get(pde_conf, :omega, 0.5)))
    elseif haskey(pde_conf, :tau)
        ConstantRelaxationTime(T(pde_conf[:tau]))
    else
        error("BGK: the PDE namespace needs a relaxation time (:tau, or :viscosity with :omega).")
    end
    nodes, weights = velocity_grid(get(pde_conf, :velocity_grid, :uniform)::Symbol, pde_conf[:Nv]::Int, pde_conf[:vmax], T)
    gamma = get(pde_conf, :gamma, (D + 2) / D)
    center = get(pde_conf, :vcenter, ntuple(_ -> 0.0, D))
    return BGKEquation(Val(D), nodes, weights, gamma, relaxation; center = center)
end
