# Time Integration & Execution Loop

The time integration module drives the physical evolution of the simulation. It handles dynamic timestep calculation, Runge-Kutta/IMEX stage execution, progress logging, and intermediate state serialization.

## The Core Orchestrator (`solve_equation`)

The `solve_equation` function acts as the primary simulation orchestrator governing the main time-stepping loop. 

*   **Initialization:** It accepts the selected Runge-Kutta or IMEX time integrator (`timestepper`), the physical equation system (`eq`), the active mesh-free domain (`pg`), the final simulation time (`tmax`), and a baseline time step (`dt`).
*   **Dynamic Time-Stepping:** If the `is_cfl` flag is true, the `dt` parameter is treated as a CFL number. The actual physical time step is then dynamically computed at each iteration using the current properties of the grid and the divergence interpolator.
*   **Progress Tracking:** The loop supports visual progress tracking (toggled via `show_progress`), which logs the simulation's advancement and dynamically calculates an ETA formatted as HH:MM:SS based on the `progress_interval`.
*   **Parallelization:** The orchestrator logs whether it is utilizing `@threads` or `@batch` for parallel execution during the integration steps.
*   **Outputs:** Upon completion of the time-stepping loop, the function returns a tuple containing the position history, state history, time history, total simulation steps (`k_step`), and the total elapsed wall time.

## State Archiving (`saveData!`)

To generate actionable output without overwhelming system memory, the framework periodically extracts the simulation state across a user-defined number of `snapshots`. This is handled exclusively by the `saveData!` function.

*   **Time Logging:** It records the exact current simulation time into the `ts_storage` array for the specific snapshot index.
*   **Memory Allocation:** The function allocates new continuous state and position vectors for the targeted snapshot.
*   **Ghost Particle Filtering:** If the `remove_ghosts` flag is set to true, the function actively filters out boundary particles and saves only the active domain core. 
*   **Unfiltered Archiving:** If `remove_ghosts` is false, it copies the corresponding views for the entire particle grid into the storage arrays.

## Documentation

```@docs
solve_equation
saveData!
```