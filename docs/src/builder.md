# The Builder Architecture

The HyperCloud framework utilizes a strictly modular factory pattern to construct simulation environments. Instead of manually instantiating nested numerical structs, users provide flat configuration dictionaries. Central orchestrators (such as `run_direct_simulation` or `run_kinetic_simulation`) segregate these dictionaries into isolated namespaces and sequentially invoke a series of `build_*` functions to assemble the computational pipeline.

## The Context Dictionary

During the assembly phase, the orchestrator maintains a shared `context::Dict{Symbol, Any}`. This dictionary acts as a dependency injection container. As each component is built, it is cached in the context (e.g., `context[:Equation]`, `context[:Domain]`, `context[:Grid]`), making it instantly accessible to downstream builders. It also stores global execution parameters, such as the spatial dimension (`context[:D]`), the precision type (`context[:Type]`), and critical spatial metrics (`context[:max_dx]`).

## Overloading and Custom Logic

The builder system is explicitly designed for user extensibility via Julia's multiple dispatch and value types (`Val`). Every major component builder extracts a specific key (usually `:name` or `:domain`) and dispatches to a targeted function signature.

If you want to inject custom physics, schemes, or boundary conditions, you do not need to modify the framework's source code. You only need to fulfill the input/output contract by overloading the generic fallback method in your own script and ensuring your struct subtypes the correct abstract interface. 

For example, to implement a custom spatial scheme, you would define your struct as a subtype of `DivergenceInterpolator`, and then overload `build_scheme`:

```julia
# 1. Define your custom scheme fulfilling the API contract
struct MyCustomScheme{T, D} <: DivergenceInterpolator{D}
    my_param::Float64
end

# 2. Overload the builder method in your local scope
function HyperCloud.build_scheme(::Val{:my_custom_scheme}, conf::Dict, context::Dict)
    T = context[:Type]::DataType
    D = context[:D]::Int
    
    # Extract custom parameters from the configuration dictionary
    my_param = conf[:my_param]
    
    # Return the instantiated object
    return MyCustomScheme{T, D}(my_param)
end
```
By passing `:name => :my_custom_scheme` in your `Scheme` configuration namespace, the framework will now dynamically route to your implementation. As long as the returned object adheres to the expected abstract type interfaces (and implements the required mathematical evaluation logic expected by the time stepper), the central orchestrator will seamlessly integrate it into the main solver loop.

## API Documentation

```@docs
build_equation
build_domain
build_particle_grid
build_ic
build_kinetic_system
build_scheme
build_flux
build_limiter
build_mood
build_mood_criterion
build_mood_strategy
build_tableau
build_timestepper
build_weights
```