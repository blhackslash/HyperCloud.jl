# Configuration Reference

A simulation is described by one flat `Dict{Symbol, Any}`. Keys are grouped into namespaces
by prefix: `:Grid_Ns` belongs to the `Grid` namespace and reaches the builder as `:Ns`
(see `extract_namespace`). Keys without a prefix are global.

A method dictionary (for example from [`create_method`](@ref)) is merged on top of the shared
parameters, so any key below may appear in either place.

## Global keys

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:sim_func_name` | — | Pipeline to run: `:run_direct_simulation` or `:run_kinetic_simulation`. Set automatically by `create_method`. |
| `:tmax` | — | Final time. `:Time_tmax` overrides it. |
| `:snapshots` | — | Number of output intervals; `snapshots + 1` states are stored, including `t = 0`. |
| `:remove_ghosts` | `true` | Drop boundary/ghost particles from the stored snapshots. |
| `:T` | `Float64` | Floating-point type of the simulation. |
| `:real_type` | `Float64` | Floating-point type of the exact reference (`analytical_solution`). |
| `:save_relax` | `false` | Kinetic runs only: store the raw kinetic variables instead of the macroscopic state. |

## `PDE_` — equation

| Key | Used by | Meaning |
|:----|:--------|:--------|
| `:name` | all | `:linear`, `:burgers` or `:euler`. Required. |
| `:D` | `:burgers`, `:euler` | Spatial dimension. Linear advection infers `D` from `:velocities`. |
| `:velocities` | `:linear` | One entry per dimension describing ``A_d``; see [Linear Advection](PDEs/advection.md). |
| `:gamma` | `:euler` | Ratio of specific heats. |

The Euler representation is read from the same namespace by `parse_representation`
(default: conservative).

## `Grid_` — domain and particles

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:domain` | — | `:rectangular`, `:spherical`, or a pre-built `GeometricDomain`. Required. |
| `:mins`, `:maxs` | — | Bounding box (tuples of length `D`). A spherical domain uses their midpoint as center and half the first extent as radius. |
| `:periodic` | `false` | `Bool` or `NTuple{D, Bool}`. |
| `:bc` | none | Boundary conditions: `Dict(tag => bc)` where `bc` is an `AbstractBoundaryCondition`, `:outflow`, `:fixed_dirichlet`, or the name of a zero-argument constructor (looked up in `Main`, then in HyperCloud). A string is parsed as Julia code. |
| `:Ns` | — | Particles per direction (tuple). The nominal spacing is `(maxs - mins) ./ Ns`. |
| `:randomness_factor` | — | Jitter per direction as a fraction of the spacing (tuple). |
| `:seed` | — | RNG seed for the jitter. Ignored (fixed to 42) when all jitter is zero. |

!!! warning "The key is `:Grid_bc`"
    Boundary conditions are read from `:Grid_bc`. A key such as `:Grid_bc_map` is silently ignored.

## `Weight_` — MLS weight function

| Key | Meaning |
|:----|:--------|
| `:name` | `:exponential` or `:inverse`. Required. |
| `:range` | Support radius in multiples of the largest grid spacing. Required. High reconstruction orders need a larger range (the examples use `5.5`). |
| `:alpha` | Shape parameter of the weight. Required. |

## `IC_` — initial condition

`:name` selects the initial condition; every other key is converted to the simulation type `T`.

| `:name` | Type | Keys |
|:--------|:-----|:-----|
| `:gauss` | [`Gauss`](@ref) | `:a` (amplitude, length `M`), `:b` (center, length `D`), `:width` |
| `:box` | [`Box`](@ref) | `:u_bg`, `:u_box` (length `M`), `:mins`, `:maxs` (length `D`) |
| `:sine` | [`Sine`](@ref) | `:a`, `:c_offset` (length `M`), `:period` (length `D`) |
| `:riemann` | [`Riemann`](@ref) | `:uL`, `:uR` (length `M`), `:p0` (length `D`), `:n` (optional, default `+x`, normalized) |
| `:s_riemann` | [`SRiemann`](@ref) | as `:riemann`, plus `:width` |
| `:q_riemann` | [`QRiemann`](@ref) | `:states` (tuple of `2^D` states, binary ordering), `:p0` |

States are given in the variables the equation evolves. For Euler these are the conservative
variables ``(\rho, \rho\vec u, E)``.

## `Scheme_` — spatial discretization

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:name` | — | `:MUSCL`, `:Upwind`, `:Central` or `:WENO`. Required. |
| `:order` | — | Formal order of the reconstruction. `:Upwind` only supports `1`. |
| `:MLS_order` | `0` | Order of the MLS divergence operator. |
| `:upwind_alg_nd` | — | `:Upwind` only: multi-dimensional upwind algorithm, e.g. `:Classic`. |

## `Flux_`, `Limiter_`, `MOOD_`

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:Flux_name` | — | `:Rusanov` or `:Upwind`. Required. |
| `:Limiter_name` | `:none` | `:none`, `:minmod`, `:superbee`, `:VK` (Venkatakrishnan), `:BJ` (Barth–Jespersen). |
| `:Limiter_mode` | — | Required by every active limiter, e.g. `:hard`. |
| `:MOOD_criterion` | `:none` | `:none`, `:only`, `:U1`, `:U2`. |
| `:MOOD_strategy` | `:none` | `:EPD0`, `:SEPD0`, `:EPD1`, `:EPD2`. Required when a criterion is set. |
| `:MOOD_delta_relax` | — | Required by `:U1`/`:U2`. Multiplied by the cell volume `prod(dx)`. |

## `Time_` — time stepping

| Key | Meaning |
|:----|:--------|
| `:stepper` | Explicit: `:Euler`, `:RK2` (Ralston), `:RK3` (SSP), `:RK4` (classical). IMEX: `:IMEXEuler`, `:ARS222`, `:ARS233`, `:PRSSP3`, `:SSP332`. |
| `:CFL` | CFL number; the step is recomputed every iteration. |
| `:dt` | Fixed time step, used when `:CFL` is absent. |
| `:tmax` | Optional override of the global `:tmax`. |

IMEX steppers treat the source terms implicitly and are used with the kinetic pipeline.

## `Kinetic_` — relaxation system

Only read by `run_kinetic_simulation`.

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:velocities` | — | Kinetic speeds, one `NK`-tuple per dimension. Use [`generate_relax_params`](@ref). |
| `:indices` | — | `Vector{Int}` of edges mapping kinetic to macroscopic variables (see [`Kin2Macro`](@ref)). |
| `:epsilon` | — | Relaxation time ``\epsilon``. |
| `:interior_factor` | `D` | Scaling of the inverse kinetic speeds in the Maxwellian. |
| `:save_relax` | `false` | Same as the global key. |

## Complete example

```julia
using HyperCloud, StaticArrays

vels, idx = generate_relax_params(2.0, 1, 1)

params = Dict{Symbol, Any}(
    :sim_func_name => :run_kinetic_simulation,
    :tmax => 0.5, :snapshots => 10,

    :PDE_name => :burgers, :PDE_D => 1,

    :Grid_domain => :rectangular, :Grid_mins => (-1.0,), :Grid_maxs => (1.0,),
    :Grid_periodic => true, :Grid_Ns => (200,),
    :Grid_randomness_factor => (0.2,), :Grid_seed => 10,

    :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 5.5,

    :IC_name => :gauss, :IC_a => (1.0,), :IC_b => (0.0,), :IC_width => 0.5,

    :Scheme_name => :MUSCL, :Scheme_order => 2,
    :Flux_name => :Rusanov, :Limiter_name => :none,
    :MOOD_criterion => :U2, :MOOD_strategy => :EPD1, :MOOD_delta_relax => 1e-4,

    :Time_stepper => :ARS222, :Time_CFL => 0.2,

    :Kinetic_velocities => vels, :Kinetic_indices => idx, :Kinetic_epsilon => 1e-6,
)

sim = run_kinetic_simulation(params)   # `nothing` if anything failed (the error is logged)
exact = analytical_solution(params)
exact(SVector(0.3, 0.5))               # exact state at x = 0.3, t = 0.5
```
