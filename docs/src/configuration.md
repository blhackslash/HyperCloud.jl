# Configuration Reference

A simulation is described by one flat `Dict{Symbol, Any}`. Keys are grouped into namespaces
by prefix: `:Grid_Ns` belongs to the `Grid` namespace and reaches the builder as `:Ns`
(see [`extract_namespace`](@ref)). Keys without a prefix are global.

Most namespaces configure numerical components of HyperCloudCore; their builders and defaults
are documented in HyperCloudCore's [Config-Driven Builders](@extref HyperCloudCore Config-Driven-Builders).

A method dictionary (for example from [`create_method`](@ref)) is merged on top of the shared
parameters, so any key below may appear in either place.

## Global keys

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:sim_func_name` | — | Pipeline to run: `:run_direct_simulation`, `:run_kinetic_simulation` or `:run_bgk_simulation`. Set automatically by `create_method`; BGK configurations (`:PDE_name => :bgk`) run [`run_bgk_simulation`](@ref) from either of the other two. |
| `:tmax` | — | Final time. `:Time_tmax` overrides it. |
| `:snapshots` | — | Number of output intervals; `snapshots + 1` states are stored, including `t = 0`. |
| `:remove_ghosts` | `true` | Drop boundary/ghost particles from the stored snapshots. |
| `:T` | `Float64` | Floating-point type of the simulation. |
| `:real_type` | `Float64` | Floating-point type of the exact reference (`analytical_solution`). |
| `:save_relax` | `false` | Kinetic and BGK runs only: store the raw kinetic variables (distributions) instead of the macroscopic state. |

## `PDE_` — equation

| Key | Used by | Meaning |
|:----|:--------|:--------|
| `:name` | all | `:linear`, `:burgers`, `:euler` or `:bgk`. Required. |
| `:D` | `:burgers`, `:euler`, `:bgk` | Spatial dimension. Linear advection infers `D` from `:velocities`. |
| `:velocities` | `:linear` | One entry per dimension describing ``A_d``; see [Linear Advection](PDEs/advection.md). |
| `:gamma` | `:euler`, `:bgk` | Ratio of specific heats (BGK default `(D + 2) / D`; it fixes the internal degrees of freedom, see [BGK](PDEs/bgk.md)). |
| `:Nv` | `:bgk` | Velocity nodes per direction (`Nv^D` discrete velocities). |
| `:vmax` | `:bgk` | Half-width of the velocity grid (outermost node for Gauss–Hermite). |
| `:velocity_grid` | `:bgk` | `:uniform` (default) or `:gauss_hermite`. |
| `:vcenter` | `:bgk` | Centre of the velocity grid per direction (default `0`). |
| `:tau` | `:bgk` | Constant relaxation time; `Inf` is collisionless. |
| `:viscosity`, `:omega` | `:bgk` | Instead of `:tau`: `τ = μ / p` with `μ = viscosity * θ^omega` (`omega` default `0.5`). |

The Euler representation is read from the same namespace by `parse_representation`:
`:representation` is `:conservative` (default), `:primitive` or `:lagrangian`, the latter two with
an integration `:path` (default `:mapped`).

## `Grid_` — domain and particles

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:domain` | — | `:rectangular`, `:spherical`, or a pre-built `GeometricDomain`. Required. |
| `:mins`, `:maxs` | — | Bounding box (tuples of length `D`). A spherical domain uses their midpoint as center and half the first extent as radius. |
| `:periodic` | `false` | `Bool` or `NTuple{D, Bool}`. |
| `:bc` | none | Boundary conditions: `Dict(tag => bc)` where `bc` is an `AbstractBoundaryCondition`, `:outflow`, `:fixed_dirichlet`, or the name of a zero-argument constructor (looked up in HyperCloudCore, then in `Main`; see [`parse_bc`](@extref HyperCloudCore HyperCloudCore.parse_bc)). A string is parsed as Julia code. |
| `:Ns` | — | Particles per direction (tuple). The nominal spacing is `(maxs - mins) ./ Ns`. |
| `:randomness_factor` | — | Jitter per direction as a fraction of the spacing (tuple). |
| `:seed` | `42` | RNG seed for the jitter. |
| `:mover` | `:none` | Grid motion: `:none`, `:lagrangian`, `:scaled_lagrangian` (factor `:mover_alpha`) or `:sinusoidal` (amplitude `:mover_amplitude`); see [Moving Grids](@extref HyperCloudCore Moving-Grids-(ALE)). Lagrangian movers need `fluid_velocity` for the equation; all built-in equations define it. |
| `:volumes` | `:shepard` | Particle volume estimator: `:shepard` (`:volumes_sweeps`) or `:mls` (`:volumes_degree`, `:volumes_samples_per_dx`); see [Particle Volumes](@extref HyperCloudCore Particle-Volumes-and-Adaptivity). Used for the initial volumes and, on a moving grid, for the re-estimation after every step and every adaptation call (HyperCloudCore ≥ 0.7.5; before, only the initial volumes). On a moving grid Shepard volumes give a resolution-independent bias in the discrete mass; `:mls` converges but costs about 3 (1-D) to 20 (2-D) time steps per re-estimation. |

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
| `:order` | — | Consistency order (expected convergence rate; MUSCL reconstructs with degree `order - 1`). `:Upwind` only supports `1`, `:WENO` needs at least `2`. |
| `:MLS_order` | `0` | Fixed order of the MLS divergence operator; `0` follows each particle's current order. |
| `:upwind_alg_nd` | — | `:Upwind` only: `:Classic` (any system), `:Tiwari` (scalar) or `:Praveen` (scalar, 2D). |

The schemes are described in HyperCloudCore: [MUSCL](@extref HyperCloudCore MUSCL-Reconstruction-and-Limiting),
[Upwind](@extref HyperCloudCore Upwind-Divergence), [Central](@extref HyperCloudCore Central-Divergence) and
[WENO](@extref HyperCloudCore Meshfree-WENO) (scalar equations only).

## `Flux_`, `Limiter_`, `MOOD_`

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:Flux_name` | — | `:Rusanov` or `:Upwind`. Required. |
| `:Limiter_name` | `:none` | `:none`, `:minmod`, `:superbee`, `:VK` (Venkatakrishnan), `:BJ` (Barth–Jespersen). |
| `:Limiter_mode` | — | Required by every active limiter, e.g. `:hard`. |
| `:MOOD_criterion` | `:none` | `:none`, `:only`, `:U1`, `:U2`. |
| `:MOOD_strategy` | `:none` | `:EPD0`, `:SEPD0`, `:EPD1`, `:EPD2`. Required when a criterion is set. |
| `:MOOD_delta_relax` | — | Required by `:U1`/`:U2`. Multiplied by the cell volume `prod(dx)`. |
| `:MOOD_plateau` | `1e-3` | `:U2` only (HyperCloudCore ≥ 0.7.6): a component that varies across the stencil by less than this fraction of its global range is not checked. Keeps smooth far fields (e.g. a vortex tail) at full order; `0` turns it off. |
| `:MOOD_curvature_ratio` | `0.5` | `:U2` only (HyperCloudCore ≥ 0.7.6): Diot's curvature-magnitude test, applied among the nearest neighbours; `0` turns it off. |

See HyperCloudCore's [Numerical Flux Functions](@extref HyperCloudCore Numerical-Flux-Functions), the
[limiters](@extref HyperCloudCore Slope-limiting-(*a-priori*)) and
[MOOD](@extref HyperCloudCore Multi-Dimensional-Optimal-Order-Detection-(MOOD)) for what the options do.

## `Time_` — time stepping

| Key | Meaning |
|:----|:--------|
| `:stepper` | Explicit: `:Euler`, `:RK2` (SSP), `:RK2Ralston`, `:RK3` (SSP), `:RK4` (classical). IMEX: `:IMEXEuler`, `:ARS222`, `:ARS233`, `:PRSSP3`, `:SSP332`. |
| `:CFL` | CFL number; the step is recomputed every iteration. It is not a classical Courant number, see [what the CFL number means on a point cloud](@extref HyperCloudCore What-the-CFL-number-means-on-a-point-cloud). |
| `:dt` | Fixed time step, used when `:CFL` is absent. |
| `:tmax` | Optional override of the global `:tmax`. |
| `:reset_orders_per_stage` | `true` restores the maximum order before every stage (default `false`: orders lowered by MOOD stay lowered for the step). |
| `:full_search_per_stage` | `true` rebuilds the neighbour lists at every stage of a moving grid (default `false`). |

IMEX steppers treat the source terms implicitly and are used with the kinetic pipeline. See
HyperCloudCore's [Butcher Tableaus](@extref HyperCloudCore Butcher-Tableaus) for the properties of each
tableau.

## `Adapt_` — particle insertion and merging

Optional; mainly for Lagrangian grids. See
[Particle insertion and merging](@extref HyperCloudCore Particle-insertion-and-merging).

| Key | Default | Meaning |
|:----|:--------|:--------|
| `:enabled` | `false` | Insert and merge particles between time steps. |
| `:max_gap`, `:min_gap` | `1.6`, `0.5` | Pair distances (in units of the spacing) above which a particle is inserted and below which a pair is merged. |
| `:clearance` | `0.7` | Free radius required around an insertion point. |
| `:samples_per_dx` | `2` | Lattice resolution of the conservative remap. |
| `:exclude_reduced_order` | `true` | No events near particles whose MOOD order was lowered. |
| `:max_events` | unlimited | At most this many insertions/merges per call (HyperCloudCore ≥ 0.7.5). |

Every adaptation call is exactly conservative with respect to `Σ V_i U_i`, also with several events per call
and, on a moving grid, across the volume re-estimation that follows (HyperCloudCore ≥ 0.7.5; with 0.7.4 and
earlier, calls with several events in 2-D/3-D and the next re-estimation were not, see
[Particle insertion and merging](@extref HyperCloudCore Particle-insertion-and-merging)).

The number of particles at the end is stored in the statistic `final_particle_count`.

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
