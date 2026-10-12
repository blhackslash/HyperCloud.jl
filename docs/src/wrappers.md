# Method Wrappers

Writing method dictionaries by hand is tedious and error-prone. These helpers create them with
consistent keys and labels.

## `create_method`

```julia
label, method = create_method(:RK2, :MUSCL, 2; limiter = :VK)
# label  == "RK2MUSCL2(VK)"
# method == Dict(:Time_stepper => :RK2, :Scheme_name => :MUSCL, :Scheme_order => 2,
#                :Limiter_name => :VK, :Limiter_mode => :hard, ...)
```

Returns `label => method`, where `method` uses the namespaced keys of the
[Configuration Reference](configuration.md).

- **Pipeline.** Unless `sim_func_name` is given, IMEX steppers (`:IMEXEuler`, `:ARS222`,
  `:ARS233`, `:PRSSP3`, `:SSP332`) select `:run_kinetic_simulation` and all others
  `:run_direct_simulation`.
- **Label.** `"<stepper><scheme><order>"`, followed in parentheses by any non-default flux,
  the limiter, and for MOOD the strategy (plus the criterion if it is not `:U2`). Examples:
  `"RK2MUSCL2"`, `"RK2MUSCL2(VK)"`, `"RK2MUSCL2(EPD1)"`, `"RK2MUSCL2(EPD2,U1)"`,
  `"RK4MUSCL3(Upwind,minmod)"`. `custom_label` replaces it.
- **Extra keywords** are stored as additional keys, e.g. `Scheme_upwind_alg_nd = :Classic` or
  `MOOD_delta_relax = 1e-4`. String values are converted to `Symbol`s.
- **`ignore`** lists parameters the method should ignore. The built-in direct methods pass the
  kinetic keys (`ignore_relax`).

The `:Upwind` scheme only exists at order 1.

## `get_master_method_dict`

Returns all built-in methods as `Dict(label => method)`:

- first-order baselines: `EulerUpwind1`, `RK2Upwind1`, `IMEXEulerUpwind1`, `ARS222Upwind1`,
  `ARS233Upwind1`, and `RK2WENO2`;
- direct second order: `RK2MUSCL2` without limiter, with `VK`, `BJ` or `minmod`, and with MOOD
  (`EPD0`, `EPD1`, `EPD2`, `SEPD0`, `U1`), plus `RK2MUSCL1`;
- kinetic second order: `ARS222MUSCL2` without limiter, with `VK` or `minmod`, and with MOOD
  (`EPD1`, `EPD2`, `SEPD0`, `U1`);
- high order: `RK4MUSCL2`–`RK4MUSCL5` without limiter, with `VK`, `minmod`, or MOOD (`EPD1`,
  `EPD2`), `RK4MUSCL4(EPD1,only)`, and kinetic `ARS233MUSCL2`–`ARS233MUSCL5` with MOOD.

Select methods by label in a `SimulationConfig`, e.g. `["RK2MUSCL2", "ARS222MUSCL2"]`.

## `generate_relax_params`

```julia
vels, idx = generate_relax_params(2.0, 2, 1)
# vels == ((2.0, -2.0, 0.0, 0.0), (0.0, 0.0, 2.0, -2.0))
# idx  == [1, 5]
```

Builds the standard kinetic velocity set: ``2D`` kinetic variables per macroscopic variable,
with speeds ``\pm v`` along each axis. It returns the speeds in the format expected by
`:Kinetic_velocities` and the edges for `:Kinetic_indices`. See
[Kinetic Relaxation](SourceTerms/kinetic.md).

## API

```@docs
create_method
get_master_method_dict
generate_relax_params
```
