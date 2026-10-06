# Method Configuration Wrappers

To prevent users from manually assembling highly nested numerical method dictionaries, HyperCloud provides programmatic factory functions. These streamline the creation of complex simulation pipelines and ensure standardized configurations across experiments.

## Core Functions

### `create_method`

Automatically generates a fully formed method dictionary using the flattened namespace API.

*   **Automatic Pipeline Routing:** If `sim_func_name` is not explicitly provided, it infers the pipeline based on the `timestepper`. If an IMEX scheme (e.g., `:ARS222`, `:ARS233`) is detected, it automatically routes to `:run_kinetic_simulation`; otherwise, it defaults to `:run_direct_simulation`.
*   **Dynamic Labeling:** Generates a standardized display string (e.g., `"RK2MUSCL2(EPD1)"`) encompassing the time stepper, gradient scheme, order, flux, limiter, and MOOD strategy, ensuring clean plot legends.

### `get_master_method_dict`

Returns a comprehensive `Dict` containing the optimal parameterizations for every standard numerical scheme natively supported by the framework. This includes:

*   Stable explicit baselines (e.g., `RK2Upwind1`).
*   Direct 2nd Order Total Variation Diminishing (TVD) schemes (e.g., `RK2MUSCL2` with `VK`, `minmod`, or `BJ` limiters).
*   High-order MOOD strategies (e.g., `RK4MUSCL4` with `U2` criteria and `EPD1`/`EPD2` strategies).
*   Stiffly-accurate kinetic relaxation formulations using IMEX integration (e.g., `ARS222MUSCL2`).

### `generate_relax_params`

A convenience builder specifically for kinetic systems. It takes a base velocity magnitude (`v_mag`), the spatial dimension (`dim`), and the macroscopic system size (`NM`) to generate symmetric kinetic wave speeds and uniform partition indices (`relax_indices`) that perfectly align with the `RelaxationSourceTerm` initialization.

## Documentation

```@docs
create_method
get_master_method_dict
generate_relax_params
```