# Simulation Execution Engine

The simulation functions serve as the central orchestrators of the HyperCloud framework. They are responsible for translating flat configuration dictionaries into fully assembled numerical contexts, managing the initialization phase, driving the main time-stepping loop, and securely packaging the final output data.

## General Orchestration Pipeline

All simulation runners in HyperCloud follow a strict, modular assembly pipeline based on a unified execution context:

1.  **Namespace Extraction:** The runner accepts a global `ParamDict` and uses `extract_namespace` to segregate parameters into isolated component dictionaries based on prefixes (e.g., separating `PDE_...`, `Grid_...`, `Time_...`, and `Scheme_...` keys).
2.  **Context Assembly:** A central execution `context` dictionary is populated sequentially. The mathematical equations are built first, yielding the spatial dimensions $D$ and system size $M$.
3.  **Geometry & Numerics:** Using the extracted namespaces, the framework builds the continuous domain, constructs the particle grid, and instantiates the spatial numerical schemes (Flux, Limiter, MOOD, and baseline Scheme).
4.  **Time Stepping & Execution:** The temporal tableau is assembled, and the time-step size `dt` is resolved (either directly or via a `CFL` condition). The main solver loop (`solve_equation`) executes until `tmax`, tracking runtime statistics and filtering out ghost particles.
5.  **Data Packaging:** The resulting spacetime vectors are packaged into a standardized `AbstractSimData` object, appending custom statistics like total runtime and step count.

---

## Direct Simulation (`run_direct_simulation`)

The `run_direct_simulation` function is the standard runner for explicit macroscopic PDEs (such as the Euler or Burgers' equations). 

*   **Source Terms:** By default, it configures the integration environment with empty tuples for both `ExplicitSources` and `ImplicitSources`.
*   **Initialization:** It applies the user-defined `InitialCondition` directly to the macroscopic physical equation across the particle grid.
*   **Data Pipeline:** The raw conservative state vectors computed by the solver are directly packaged into the final data structure without modification.

---

## Kinetic Simulation (`run_kinetic_simulation`)

The `run_kinetic_simulation` function manages the execution of highly coupled kinetic relaxation systems. It introduces a two-tier equation architecture to handle the bridging between macroscopic physics and kinetic transport.

*   **Two-Tier Assembly:** It first builds the base macroscopic equation (`eq_macro`) before extracting an additional `Kinetic` namespace to construct the kinetic system (`eq_kin`) and its associated `source_term`.
*   **Dimensionality Overwrite:** To seamlessly reuse the standard spatial schemes, the execution context is dynamically overwritten to operate in the higher $NK$-dimensional kinetic space, and the relaxation source term is injected into `ImplicitSources`.
*   **Kinetic Initialization:** It delegates initialization to the specialized kinetic overload of `set_initial_conditions!`, which projects the macroscopic initial condition onto the kinetic Maxwellian state.
*   **Macroscopic Collapse:** Because users typically want to analyze the macroscopic physics, the runner automatically collapses the $NK$-dimensional kinetic data back into the $NM$-dimensional macroscopic state using the source term's `km` (Kin2Macro) mapping before returning the data. This behavior can be bypassed to save the raw kinetic variables by setting `save_relax` to `true`.

---

## API Documentation

```@docs
run_direct_simulation
run_kinetic_simulation
extract_namespace
```