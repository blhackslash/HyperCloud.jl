# Simulation Engine

The runners turn a flat parameter dictionary into a complete solver, run it, and package the
result: [`run_direct_simulation`](@ref) for the PDE itself, [`run_kinetic_simulation`](@ref) for its
kinetic relaxation approximation, and [`run_bgk_simulation`](@ref) for the discrete-velocity
[BGK model](PDEs/bgk.md).

## Pipeline

Both runners follow the same steps, sharing one `context` dictionary (see
[Builder API](builder.md); the numerical builders are HyperCloudCore's
[config-driven builders](@extref HyperCloudCore Config-Driven-Builders)):

1. **Split the parameters** into namespaces with `extract_namespace` (`PDE`, `Grid`, `Weight`,
   `Flux`, `Limiter`, `MOOD`, `Scheme`, `IC`, `Time`, `Adapt`, and `Kinetic` for the kinetic runner).
2. **Build the equation.** It determines the dimension `D` and system size `M`.
3. **Build the geometry and particle grid**, including the MLS weight function.
4. **Build the spatial discretization:** flux, limiter, MOOD and scheme.
5. **Build the time stepper**, HyperCloudCore's
   [`GeneralRKTimeStepper`](@extref HyperCloudCore HyperCloudCore.GeneralRKTimeStepper), from the
   Butcher tableau (explicit RK or IMEX), and optionally the particle adaptivity (`Adapt` namespace).
6. **Apply the initial condition.**
7. **Integrate** with [`solve_equation`](@ref) up to `tmax`, using either a fixed `dt` or a CFL
   number.
8. **Package** the snapshots into an `AbstractSimData` and add the statistics `Runtime`
   (wall-clock seconds), `k_step` (number of steps) and `final_particle_count`, plus
   `reserved_particles` and `reserved_bytes` with a [memory limit](@ref "Memory limit").

The simulation precision is `params[:T]` (default `Float64`).

The returned object holds one entry per snapshot: `sim.x[k]` are the positions, `sim.u[k]`
the states and `sim.t[k]` the time of snapshot `k`.

!!! warning "Errors are caught"
    Both runners catch every exception, log it with `@error` (including the backtrace), and
    return `nothing`. Always check the result, e.g. `@assert sim !== nothing`, when calling them
    directly.

## Running on a GPU

The runners place a run where HyperCloudCore's execution backend runs ([`storage_for`](@ref);
see [Execution backends](@extref HyperCloudCore Execution-backends-and-the-GPU-roadmap)).
The CPU backends ([`PolyesterBackend()`](@extref HyperCloudCore HyperCloudCore.PolyesterBackend), the
default, [`ThreadsBackend()`](@extref HyperCloudCore HyperCloudCore.ThreadsBackend) and
`KernelAbstractions.CPU()`) keep the run in host memory. With a GPU backend, the grid, the scheme
and the time stepper live on the device (`context[:Storage]`). No parameter changes, so a GPU run
has the same parameter hash as its CPU counterpart:

```julia
using HyperCloud, AMDGPU          # or CUDA / Metal / oneAPI
set_backend!(ROCBackend())
sim = run_direct_simulation(params)
set_backend!(PolyesterBackend())  # back to the CPU

with_backend(() -> run_direct_simulation(params), ROCBackend())   # the same, scoped
```

- The grid is constructed on the host and then moved to the device; the time steps (neighbour
  search, MOOD, boundary conditions, CFL estimate) run there.
- Initial conditions are evaluated on the host and copied over once; each snapshot copies the
  active particles back. The returned `sim` therefore holds plain host arrays.
- Particle adaptivity works; a step that actually inserts or merges particles copies the grid to the
  host and back.
- All schemes, including `:WENO`, and all grid movers run on the device (since HyperCloudCore 0.6.1).
  What still runs on the host is listed under
  [Not yet on the device](@extref HyperCloudCore Not-yet-on-the-device).
- Custom equations, boundary conditions and source terms run inside the time step and must
  follow the [rules for kernel bodies](@extref HyperCloudCore Rules-for-kernel-bodies) (bits types,
  no allocation, no dynamic dispatch). Initial conditions only run on the host.

### Memory limit

By default the buffers of a run grow on demand; growing a buffer briefly holds the old and the new
array at once, which is where a run on a nearly full device fails. `set_memory_limit!(bytes)`
instead sizes every buffer of the following runs once, to the largest particle capacity under the
limit (HyperCloudCore's [`reserve!`](@extref HyperCloudCore HyperCloudCore.reserve!)), and a run that
outgrows it stops with a [`CapacityError`](@extref HyperCloudCore HyperCloudCore.CapacityError):

```julia
set_memory_limit!(12 * 2^30; slots_per_particle = 40)   # 12 GiB; room for neighbourhoods to grow
sim = run_direct_simulation(params)                       # stats: reserved_particles, reserved_bytes
set_memory_limit!(nothing)
```

Like the backend it is a session setting, not a parameter, so it does not change the run hash.
`slots_per_particle` defaults to what the initial neighbour lists need; raise it for moving or adaptive
grids.

!!! warning "Device reservations (AMDGPU 2.9.1, ROCm 7.2)"
    A reservation on a device can zero grid buffers when earlier device runs happened in the same
    Julia session (see HyperCloudCore's [memory budget](@extref HyperCloudCore Memory-budget) notes).
    The runners compare the grid before and after the reservation and refuse to run when it changed.
    Until this is fixed upstream, start reserved GPU runs in a fresh session.

## Direct simulation

[`run_direct_simulation`](@ref) evolves the PDE given in the `PDE` namespace.

- No source terms are attached.
- The initial condition is evaluated at every particle.
- The stored states are the variables the equation evolves (conservative variables for Euler).

## Kinetic simulation

[`run_kinetic_simulation`](@ref) replaces the macroscopic system with ``M`` variables by a
linear kinetic system with ``N_K`` variables plus a stiff relaxation term (see
[Kinetic Relaxation](SourceTerms/kinetic.md)).

- **Two equations.** The macroscopic equation is built first; [`build_kinetic_system`](@ref)
  then builds the kinetic `LinearAdvection` system and the [`RelaxationSourceTerm`](@ref) from
  the `Kinetic` namespace.
- **Kinetic context.** The context is switched to the kinetic equation and `M = N_K`, so the
  same schemes run on the kinetic variables. The relaxation term is registered as an implicit
  source, which requires an IMEX stepper (`:IMEXEuler`, `:ARS222`, `:ARS233`, `:PRSSP3`, `:SSP332`).
- **Initialization.** The macroscopic initial condition is mapped to the kinetic equilibrium
  (Maxwellian) at every particle.
- **Output.** Snapshots are summed back to the macroscopic variables with the source term's
  `Kin2Macro` map. Set `:save_relax => true` to keep the raw kinetic variables instead.

## BGK simulation

[`run_bgk_simulation`](@ref) runs the [BGK equation](PDEs/bgk.md) (`:PDE_name => :bgk`). The direct
and kinetic runners hand BGK configurations over to it, so the method dictionaries work unchanged.

- **Macroscopic initial data.** The initial condition is built with `M = D + 2` (the Euler variables
  of the same ``\gamma``) and every particle starts in the discrete Maxwellian of its state.
- **Collision term.** [`BGKRelaxation`](@ref) (implicit) for an IMEX tableau,
  [`ExplicitBGKRelaxation`](@ref) for an explicit one, none for `:PDE_tau => Inf`.
- **Output.** Snapshots hold the moments ``(\rho, \rho\vec u, E)``; `:save_relax => true` keeps the
  distributions.

## API

```@docs
run_direct_simulation
run_kinetic_simulation
extract_namespace
storage_for
with_backend
set_memory_limit!
get_memory_limit
```
