# Builder API

The runners never construct solver components directly. They split the parameter dictionary
into namespaces with [`extract_namespace`](@ref) and call one `build_*` function per component.
Every builder dispatches on a `Val` of a configuration value, so new components are added by
defining a method, without changing HyperCloud.

The builders come from two packages:

- **Physics, in HyperCloud:** [`build_equation`](@ref), [`build_ic`](@ref) and
  [`build_kinetic_system`](@ref) (documented below).
- **Numerics, in HyperCloudCore:** domain, particle grid, weights, flux, limiter, MOOD, scheme,
  tableau, time stepper and adaptivity. They are re-exported by HyperCloud and documented in
  HyperCloudCore's [Config-Driven Builders](@extref HyperCloudCore Config-Driven-Builders),
  together with the [context keys](@extref HyperCloudCore Context-keys) they read and write.

## Build order and context

The builders share a `context::Dict{Symbol, Any}`. Each one reads what earlier builders stored
and adds its own result. [`run_direct_simulation`](@ref) runs them in this order:

| Step | Builder | Reads from `context` | Stores |
|:-----|:--------|:---------------------|:-------|
| 0 | runner | — | `:T`, `:Storage` ([`storage_for`](@ref)), `:MemoryLimit` (only with [`set_memory_limit!`](@ref)) |
| 1 | [`build_equation`](@ref) | `:T` | runner stores `:Equation`, `:D`, `:M`, `:ExplicitSources`, `:ImplicitSources` |
| 2 | [`build_domain`](@extref HyperCloudCore HyperCloudCore.build_domain) | `:T`, `:D` | runner stores `:Domain`, `:WeightConf` |
| 3 | [`build_particle_grid`](@extref HyperCloudCore HyperCloudCore.build_particle_grid) | `:T`, `:D`, `:M`, `:Domain`, `:WeightConf`, `:Storage` | `:max_dx`, `:Grid` |
| 3a | [`build_weights`](@extref HyperCloudCore HyperCloudCore.build_weights), [`build_mover`](@extref HyperCloudCore HyperCloudCore.build_mover), [`build_volume_estimator`](@extref HyperCloudCore HyperCloudCore.build_volume_estimator) (called by step 3) | `:T`, `:max_dx`, `:Domain` | `:interp_range` |
| 3b | [`build_adaptivity`](@extref HyperCloudCore HyperCloudCore.build_adaptivity) | `:T` | — |
| 4 | [`build_flux`](@extref HyperCloudCore HyperCloudCore.build_flux), [`build_limiter`](@extref HyperCloudCore HyperCloudCore.build_limiter), [`build_mood`](@extref HyperCloudCore HyperCloudCore.build_mood) | `:T`, `:Grid` (MOOD) | runner stores `:Flux`, `:Limiter`, `:MOOD` |
| 5 | [`build_scheme`](@extref HyperCloudCore HyperCloudCore.build_scheme) | `:T`, `:D`, `:M`, `:Flux`, `:Limiter`, `:Storage` | runner stores `:Scheme` |
| 6 | [`build_tableau`](@extref HyperCloudCore HyperCloudCore.build_tableau), [`build_timestepper`](@extref HyperCloudCore HyperCloudCore.build_timestepper) | everything above | runner stores `:Tableau`; `:Memory` (with a memory limit) |
| 7 | [`build_ic`](@ref) | `:T`, `:D`, `:M` (macroscopic) | — |

The particle grid is built on a CPU backend and then moved to the device when `:Storage` is a
device array type, and the time stepper's memory reservation is checked against a fingerprint of
the grid (see [Running on a GPU](@ref)).

[`run_kinetic_simulation`](@ref) builds the initial condition right after step 1, while `:M`
still holds the macroscopic size. Then [`build_kinetic_system`](@ref) overwrites `:Equation`
and `:M` with the kinetic system and puts the relaxation term into `:ImplicitSources`, and the
remaining steps run on the kinetic variables.

`:T` is the floating-point type (`params[:T]`, default `Float64`). The runners set
`:ExplicitSources` and `:ImplicitSources` to `()` unless a source term is added.

## Extending the builders

| Builder | Package | Dispatches on | Extend with |
|:--------|:--------|:--------------|:------------|
| `build_equation` | HyperCloud | `:PDE_name` | `HyperCloud.build_equation(::Val{:x}, conf, ctx)` |
| `build_ic` | HyperCloud | `:IC_name` | `HyperCloud.build_ic(::Val{:x}, conf, ctx)` |
| `build_domain` | HyperCloudCore | `:Grid_domain` | `HyperCloudCore.build_domain(::Val{:x}, conf, ctx)` |
| `build_weights` | HyperCloudCore | `:Weight_name` | `HyperCloudCore.build_weights(::Val{:x}, conf, ctx)` |
| `build_mover` | HyperCloudCore | `:Grid_mover` | `HyperCloudCore.build_mover(::Val{:x}, conf, ctx)` |
| `build_volume_estimator` | HyperCloudCore | `:Grid_volumes` | `HyperCloudCore.build_volume_estimator(::Val{:x}, conf, ctx)` |
| `build_flux` | HyperCloudCore | `:Flux_name` | `HyperCloudCore.build_flux(::Val{:x}, conf, ctx)` |
| `build_limiter` | HyperCloudCore | `:Limiter_name` | `HyperCloudCore.build_limiter(::Val{:x}, conf, ctx)` |
| `build_mood_criterion` | HyperCloudCore | `:MOOD_criterion` | `HyperCloudCore.build_mood_criterion(::Val{:x}, conf, ctx)` |
| `build_mood_strategy` | HyperCloudCore | `:MOOD_strategy` | `HyperCloudCore.build_mood_strategy(::Val{:x}, conf, ctx)` |
| `build_scheme` | HyperCloudCore | `:Scheme_name` | `HyperCloudCore.build_scheme(::Val{:x}, conf, ctx)` |
| `build_tableau` | HyperCloudCore | `:Time_stepper` | `HyperCloudCore.build_tableau(::Val{:x}, conf, ctx)` |

Because HyperCloud re-exports HyperCloudCore, `HyperCloud.build_scheme` and
`HyperCloudCore.build_scheme` are the same function; qualifying with the owning package makes it
clear where the method goes. Unknown names raise an error naming the missing method. A complete
example for initial conditions is in [Initial Conditions](ICs.md); a new equation is described
in [PDE Architecture](PDEs/main.md).

A custom spatial scheme follows the same pattern. The returned object must be a
[`DivergenceInterpolator`](@extref HyperCloudCore HyperCloudCore.DivergenceInterpolator) that
implements HyperCloudCore's [scheme contract](@extref HyperCloudCore A-new-spatial-scheme-(DivergenceInterpolator)):

```julia
function HyperCloudCore.build_scheme(::Val{:my_scheme}, conf::Dict, context::Dict)
    T = context[:T]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    return MyScheme(T, D, M, conf[:my_param]; flux = context[:Flux])
end
```

`:Scheme_name => :my_scheme, :Scheme_my_param => 3` then selects it.

### Custom boundary conditions

`:Grid_bc` maps boundary tags to boundary conditions
([`parse_bc`](@extref HyperCloudCore HyperCloudCore.parse_bc)). Besides instances and the names
`:outflow` and `:fixed_dirichlet`, it accepts the name of any zero-argument constructor defined
in HyperCloudCore or in `Main`:

```julia
struct MyWall <: AbstractBoundaryCondition end   # plus the functor HyperCloudCore requires
params[:Grid_bc] = Dict(1 => :MyWall, 2 => :outflow)
```

The functor a boundary condition implements is described in HyperCloudCore's
[Boundary Conditions](@extref HyperCloudCore Boundary-Conditions).

## API

```@docs
build_equation
build_ic
build_kinetic_system
```
