# Kinetic Relaxation

The kinetic relaxation approach replaces a nonlinear system of ``N_M`` conservation laws by
``N_K`` linear transport equations with a stiff source:

```math
\partial_t f_k + \sum_{d=1}^D \lambda_{k,d}\,\partial_{x_d} f_k = \frac{\mathcal{M}_k(u) - f_k}{\epsilon},
\qquad u = \text{sum of the kinetic variables belonging to each macroscopic variable}.
```

As ``\epsilon \to 0`` the macroscopic variables ``u`` satisfy the original conservation law.
The transport part is a diagonal `LinearAdvection` system, so all spatial schemes work on it
unchanged. The source is solved implicitly.

## Setting it up

Use the kinetic pipeline with an IMEX time stepper and the `Kinetic_` namespace:

```julia
vels, idx = generate_relax_params(2.0, D, M)   # speeds ±2 along each axis

params[:sim_func_name]      = :run_kinetic_simulation
params[:Time_stepper]       = :ARS222
params[:Kinetic_velocities] = vels
params[:Kinetic_indices]    = idx
params[:Kinetic_epsilon]    = 1e-6
```

The kinetic speed must exceed the largest wave speed of the macroscopic system
(subcharacteristic condition). For Burgers with ``|u| \le 1``, a speed of `2.0` is safe.

## Mapping kinetic to macroscopic variables (`Kin2Macro`)

[`Kin2Macro`](@ref) is built from a list of edges. `Kin2Macro([1, 3, 6])` assigns kinetic
variables `1:2` to macroscopic variable 1 and `3:5` to variable 2.

- `km(f)` sums each group and returns the macroscopic state.
- `km(k)` returns the macroscopic index of kinetic variable `k`.

[`generate_relax_params`](@ref)`(v, D, M)` produces ``2D`` variables per macroscopic variable,
with speeds ``\pm v`` along each axis, and the matching edges.

## Equilibrium and source term

For kinetic variable ``k`` belonging to macroscopic variable ``m`` the equilibrium is

```math
\mathcal{M}_k(u) = c_m\Big(u_m + \sum_{d=1}^D F_d^m(u)\, s_{k,d}\Big),
\qquad
s_{k,d} = \begin{cases} \gamma / \lambda_{k,d} & \lambda_{k,d} \neq 0 \\ 0 & \text{otherwise,}\end{cases}
```

where ``c_m = 1 / (\text{number of kinetic variables of } m)``, ``F_d^m`` is the physical flux,
and ``\gamma`` is `interior_factor` (default ``D``). The constructor precomputes the ``s_{k,d}``
and ``1/\epsilon``.

For the velocity sets from `generate_relax_params` with the default `interior_factor`, the
equilibrium reproduces the macroscopic state and flux exactly:

```math
\sum_{k \in m} \mathcal{M}_k(u) = u_m,
\qquad
\sum_{k \in m} \lambda_{k,d}\,\mathcal{M}_k(u) = F_d^m(u).
```

`evaluate_source` returns ``(\mathcal{M}_k(u) - f_k)/\epsilon``.

## Implicit solve

In each IMEX stage with implicit coefficient ``a_{ii}``, the stage value satisfies
``Y = Y^* + a_{ii}\Delta t\,(\mathcal{M}(u) - Y)/\epsilon``. The relaxation does not change the
macroscopic variables (``\sum_{k\in m}(\mathcal{M}_k - f_k) = 0``), so ``u`` can be taken from
``Y^*``, and `implicit_solve` is exact and explicit in form:

```math
Y_k = \frac{Y^*_k + \frac{a_{ii}\Delta t}{\epsilon}\,\mathcal{M}_k(u(Y^*))}{1 + \frac{a_{ii}\Delta t}{\epsilon}} .
```

This holds for any ``\epsilon > 0``, including the stiff limit.

## Initialization

[`set_initial_conditions!`](@ref) has a kinetic method. It takes the macroscopic initial
condition and sets every particle to the equilibrium ``\mathcal{M}(u_0(x))``, so you never write
kinetic initial data by hand.

## API

```@docs
Kin2Macro
RelaxationSourceTerm
```
