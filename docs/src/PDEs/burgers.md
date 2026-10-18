# Burgers' Equation

`BurgersEquation{D, T}` implements the inviscid Burgers' equation in ``D`` dimensions:

```math
\frac{\partial u}{\partial t} + \sum_{d=1}^D \frac{\partial}{\partial x_d}\left(\frac{u^2}{2}\right) = 0 .
```

The state is a scalar (`M = 1`). Information travels with velocity ``(u, \dots, u)``, so
compressive regions steepen into shocks and expansive regions open into rarefaction fans.

In a configuration: `:PDE_name => :burgers, :PDE_D => D`.

## Implementation

- **Variables.** `prim2cons` and `cons2prim` are the identity.
- **Flux.** ``F_d(u) = u^2/2`` in every direction.
- **Wave speed.** `max_eigenvalue` returns ``|u|``; `velocity` returns the ``1 \times 1``
  Jacobian ``u``.
- **Moving grids.** `fluid_velocity` returns ``(u, \dots, u)``, the characteristic velocity, so a
  Lagrangian grid follows the characteristics.

## Exact solutions

### Riemann and smoothed Riemann (any dimension)

Along the normal coordinate ``\xi = (x - p_0)\cdot n`` the equation reduces to a 1D law with
flux ``a\,u^2/2``, where ``a = \sum_d n_d``. The characteristic speeds are ``a u_L`` and
``a u_R``:

- ``a u_L > a u_R``: a shock moving at ``a (u_L + u_R)/2`` (Rankine–Hugoniot).
- otherwise: a rarefaction fan with ``u = \xi / (a t)`` between the two speeds.

The wave type depends on the sign of ``a``, so flipping the normal and swapping the states gives
the same solution. The normal does not need to be normalized. For `SRiemann` the closure returns
the sharp Riemann solution, a good reference once the waves are wider than `width`.

### Gauss (any dimension)

The solution is given implicitly by ``u = u_0(x - u\,t\,\mathbf{1})``, which is solved with
Newton's method. This is only valid before the pulse breaks, at

```math
t^* = \frac{\mathrm{width}}{|a|}\sqrt{\frac{e}{2D}} .
```

Past ``t^*`` Newton may not converge; the closure then warns and returns its last iterate.

### Box (1D)

The top-hat case (`u_box > u_bg`) produces a rarefaction at the left edge and a shock at the
right edge. The well case (`u_box < u_bg`) produces a shock at the left edge and a rarefaction at
the right edge. Both are exact before and after the time the fan catches up with the shock,
``t_\mathrm{int} = 2(x_e - x_s)/|u_\mathrm{box} - u_\mathrm{bg}|``. After that the shock moves like
``\sqrt{t}`` relative to the background flow and its strength decays like ``t^{-1/2}``.

### Other initial conditions

There is no exact solution. The closure warns and returns the initial condition, which is
only correct at ``t = 0``.

## API

```@docs
BurgersEquation
```
