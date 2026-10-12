# Simulation Engine

The runners turn a flat parameter dictionary into a complete solver, run it, and package the
result. There are two: [`run_direct_simulation`](@ref) for the PDE itself and
[`run_kinetic_simulation`](@ref) for its kinetic relaxation approximation.

## Pipeline

Both runners follow the same steps, sharing one `context` dictionary (see
[Builder API](builder.md)):

1. **Split the parameters** into namespaces with `extract_namespace` (`PDE`, `Grid`, `Weight`,
   `Flux`, `Limiter`, `MOOD`, `Scheme`, `IC`, `Time`, and `Kinetic` for the kinetic runner).
2. **Build the equation.** It determines the dimension `D` and system size `M`.
3. **Build the geometry and particle grid**, including the MLS weight function.
4. **Build the spatial discretization:** flux, limiter, MOOD and scheme.
5. **Build the time stepper** from the Butcher tableau (explicit RK or IMEX).
6. **Apply the initial condition.**
7. **Integrate** with [`solve_equation`](@ref) up to `tmax`, using either a fixed `dt` or a CFL
   number.
8. **Package** the snapshots into an `AbstractSimData` and add the statistics `Runtime`
   (wall-clock seconds) and `k_step` (number of steps).

The simulation precision is `params[:T]` (default `Float64`).

The returned object holds one entry per snapshot: `sim.x[k]` are the positions, `sim.u[k]`
the states and `sim.t[k]` the time of snapshot `k`.

!!! warning "Errors are caught"
    Both runners catch every exception, log it with `@error` (including the backtrace), and
    return `nothing`. Always check the result, e.g. `@assert sim !== nothing`, when calling them
    directly.

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

## API

```@docs
run_direct_simulation
run_kinetic_simulation
extract_namespace
```
