# HyperCloud.jl

**HyperCloud.jl** is a meshfree framework for hyperbolic conservation laws. It solves PDEs
on particle clouds with moving-least-squares (MLS) reconstructions: first-order upwind, MUSCL
up to fifth order, and WENO. Slope limiters and a posteriori MOOD limiting keep shocks stable.
A kinetic relaxation pipeline with IMEX time stepping approximates nonlinear systems by linear
transport plus a stiff source term.

HyperCloud is part of the **PDEStudio** ecosystem:

- [HyperCloudCore.jl](@extref HyperCloudCore HyperCloudCore:-A-Meshfree-Hyperbolic-Solver)
  provides the numerical engine (particle grids, MLS schemes, MOOD, time steppers, boundary
  conditions, GPU backends) and the [config-driven builders](@extref HyperCloudCore Config-Driven-Builders)
  that HyperCloud's runners call. HyperCloud re-exports it.
- `PDEStudioCore.jl` runs batches of simulations, caches results on disk under a hash of
  their parameters, and computes error statistics.
- `PDEStudio.jl` plots the results, either locally (`GLMakie`) or in the browser from a
  headless machine (`WGLMakie` + Bonito).

## A first simulation

A simulation is a function returning a `SimulationConfig`: shared parameters, a dictionary of
numerical methods, and the methods to run.

```julia
using HyperCloud

function my_advection()
    shared = Dict{Symbol, Any}(
        :sim_func_name => :run_direct_simulation,
        :tmax => 1.0, :snapshots => 20,

        :PDE_name => :linear, :PDE_velocities => ((1.0,),),   # u_t + u_x = 0

        :Grid_domain => :rectangular, :Grid_mins => (0.0,), :Grid_maxs => (1.0,),
        :Grid_Ns => (200,), :Grid_periodic => true,
        :Grid_randomness_factor => (0.2,), :Grid_seed => 42, :Grid_mover => :none,

        :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 5.5,

        :IC_name => :gauss, :IC_a => (1.0,), :IC_b => (0.5,), :IC_width => 0.1,

        :Time_CFL => 0.2,
    )

    methods = Dict{Symbol, Dict{Symbol, Any}}(
        :analytical_solution => Dict{Symbol, Any}(),
        :RK2MUSCL2 => last(create_method(:RK2, :MUSCL, 2)),
        :RK4MUSCL4 => last(create_method(:RK4, :MUSCL, 4)),
    )

    return SimulationConfig(shared, methods, [:analytical_solution, :RK2MUSCL2, :RK4MUSCL4];
                            ref_func_name = :analytical_solution)
end

@sim my_advection     # run and cache
@plot my_advection    # open the interactive viewer (needs PDEStudio and a Makie backend)
```

Every available key is listed in the [Configuration Reference](configuration.md).

## Workflow

1. **Configure.** Write a function returning a `SimulationConfig`, either in your script or as
   `name.jl` in one of the `CONFIG_SEARCH_PATHS` folders (see [Execution Macros](macros.md)).
2. **Run.** `@sim name` builds and runs every active method and caches the results.
3. **Visualize.** `@plot name` loads the configuration into the PDEStudio GUI.
4. **Archive.** `@backup name` copies the configuration file into a timestamped `backups` folder.

## Package layout

| Folder / file | Content |
|:--------------|:--------|
| `PDEs/` | Equation definitions (linear advection, Burgers, Euler) and their exact solutions. |
| `SourceTerms/` | Source terms, currently the kinetic relaxation term. |
| `ICs/` | Initial conditions: Gaussian pulse, box, sine wave and three Riemann variants. |
| `Builder/` | The physics builders (`build_equation`, `build_kinetic_system`) and the exact references; the numerical builders come from HyperCloudCore (see [Builder API](builder.md)). |
| `device.jl` | GPU placement and memory limit of the runners ([`storage_for`](@ref), [`set_memory_limit!`](@ref)). |
| `inbuilt_methods.jl` | [`create_method`](@ref), [`get_master_method_dict`](@ref), [`generate_relax_params`](@ref). |
| `SimulationFunctions/` | The pipelines [`run_direct_simulation`](@ref) and [`run_kinetic_simulation`](@ref). |
| `time_integration.jl` | The time loop [`solve_equation`](@ref) and snapshot storage. |
| `Builder/reference.jl` | [`analytical_solution`](@ref) and the [`analytic_closure`](@ref) interface. |
| `macros.jl` | [`@sim`](@ref), [`@plot`](@ref), [`@backup`](@ref). |
| `Examples/` | Ready-to-run benchmarks: `euler_1d_sod`, `burgers_kinetic`, `advection_2d`. |
