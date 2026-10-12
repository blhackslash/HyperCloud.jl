# Builder API

The runners never construct solver components directly. They split the parameter dictionary
into namespaces and call one `build_*` function per component. Every builder dispatches on a
`Val` of a configuration value, so new components are added by defining a method, without
changing HyperCloud.

## Build order and context

The builders share a `context::Dict{Symbol, Any}`. Each one reads what earlier builders stored
and adds its own result:

| Step | Builder | Reads from `context` | Stores |
|:-----|:--------|:---------------------|:-------|
| 0 | runner | — | `:T` |
| 1 | [`build_equation`](@ref) | `:T` | `:Equation`, `:D`, `:M` |
| 1b | [`build_kinetic_system`](@ref) (kinetic only) | `:T`, `:D`, `:M`, `:Equation` | overwrites `:Equation`, `:M`; `:ImplicitSources` |
| 2 | [`build_domain`](@ref) | `:T`, `:D` | `:Domain` |
| 3 | [`build_particle_grid`](@ref) | `:T`, `:D`, `:M`, `:Domain`, `:WeightConf` | `:max_dx`, `:Grid` |
| 3b | [`build_weights`](@ref) (called by step 3) | `:T`, `:max_dx` | `:interp_range` |
| 4 | [`build_flux`](@ref), [`build_limiter`](@ref), [`build_mood`](@ref) | `:T`, `:Grid` (MOOD) | `:Flux`, `:Limiter`, `:MOOD` |
| 5 | [`build_scheme`](@ref) | `:T`, `:D`, `:M`, `:Flux`, `:Limiter` | `:Scheme` |
| 6 | [`build_tableau`](@ref), [`build_timestepper`](@ref) | everything above, `:ExplicitSources`, `:ImplicitSources` | `:Tableau` |
| 7 | [`build_ic`](@ref) | `:T`, `:D`, `:M` (macroscopic) | — |

The kinetic runner builds the initial condition right after step 1, while `:M` still holds the
macroscopic size.

`:T` is the floating-point type (`params[:T]`, default `Float64`). The runners set
`:ExplicitSources` and `:ImplicitSources` to `()` unless a source term is added.

## Extending the builders

| Builder | Dispatches on | Extend with |
|:--------|:--------------|:------------|
| `build_equation` | `:PDE_name` | `build_equation(::Val{:x}, conf, ctx)` |
| `build_domain` | `:Grid_domain` | `build_domain(::Val{:x}, conf, ctx)` |
| `build_weights` | `:Weight_name` | `build_weights(::Val{:x}, conf, ctx)` |
| `build_ic` | `:IC_name` | `build_ic(::Val{:x}, conf, ctx)` |
| `build_flux` | `:Flux_name` | `build_flux(::Val{:x}, conf, ctx)` |
| `build_limiter` | `:Limiter_name` | `build_limiter(::Val{:x}, conf, ctx)` |
| `build_mood_criterion` | `:MOOD_criterion` | `build_mood_criterion(::Val{:x}, conf, ctx)` |
| `build_mood_strategy` | `:MOOD_strategy` | `build_mood_strategy(::Val{:x}, conf, ctx)` |
| `build_scheme` | `:Scheme_name` | `build_scheme(::Val{:x}, conf, ctx)` |
| `build_tableau` | `:Time_stepper` | `build_tableau(::Val{:x}, conf, ctx)` |

Unknown names raise an error naming the missing method. A complete example for initial
conditions is in [Initial Conditions](ICs.md); a new equation is described in
[PDE Architecture](PDEs/main.md).

A custom spatial scheme follows the same pattern. The returned object must implement the
interface the time steppers expect from the built-in schemes:

```julia
function HyperCloud.build_scheme(::Val{:my_scheme}, conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    return MyScheme(T, D, M, conf[:my_param]; flux = context[:Flux])
end
```

`:Scheme_name => :my_scheme, :Scheme_my_param => 3` then selects it.

### Custom boundary conditions

`:Grid_bc` maps boundary tags to boundary conditions. Besides instances and the names
`:outflow` and `:fixed_dirichlet`, it accepts the name of any zero-argument constructor defined
in `Main` or in HyperCloud:

```julia
struct MyWall <: AbstractBoundaryCondition end   # plus the methods HyperCloudCore requires
params[:Grid_bc] = Dict(1 => :MyWall, 2 => :outflow)
```

## API

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
