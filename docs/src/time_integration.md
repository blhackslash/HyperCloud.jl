# Time Integration

[`solve_equation`](@ref) runs the time loop: it computes the step size, advances the particle
grid with the chosen Runge–Kutta or IMEX stepper, reports progress, and stores snapshots.

## Time step

- **Fixed:** with `is_cfl = false`, `dt` is the physical time step.
- **CFL:** with `is_cfl = true` (`:Time_CFL` in a configuration), `dt` is a CFL number and the
  physical step is recomputed every iteration as `dt * get_time_step(pg, eq, div_interp)`. This
  uses the current wave speeds of the equation and the grid spacing.

The last step is shortened so the run ends exactly at `tmax`.

## Snapshots

`snapshots + 1` states are stored, at the target times `range(0, tmax, length = snapshots + 1)`.

- Snapshot 1 is the initial state at `t = 0`.
- Each further snapshot is taken at the end of the first step that reaches its target time,
  so `t[k]` can lie slightly after the target.
- The final state is always stored in the last slot.

[`saveData!`](@ref) copies positions and states into new arrays, so later steps do not modify
stored snapshots. With `remove_ghosts = true` only non-boundary particles are kept. The runners
pass `:remove_ghosts` from the configuration (default `true`); `solve_equation` itself defaults
to `false`.

## Progress output

- `show_progress = true`: a progress bar.
- `show_progress = false`: an `@info` line with the percentage and an ETA (HH:MM:SS) every
  `progress_interval` seconds, which is better suited for log files.

The runners take both settings from global switches: `enable_progress_bar!(false)` and
`set_progress_interval!(seconds)`.

The loop also logs whether parallel loops use `@threads` or `@batch`.

## Return value

`solve_equation` returns `(positions, states, times, k_step, elapsed_time)`: the stored
positions, states and times per snapshot, the number of steps, and the wall-clock time of the
loop in seconds.

## API

```@docs
solve_equation
saveData!
```
