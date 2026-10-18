# Euler Equations

`EulerEquation{D, M, T, R}` implements the compressible Euler equations for an ideal gas in
``D`` dimensions, with conservative variables ``U = (\rho, \rho\vec u, E)`` and ``M = D + 2``:

```math
\begin{aligned}
\partial_t \rho + \nabla\cdot(\rho\vec u) &= 0,\\
\partial_t (\rho\vec u) + \nabla\cdot(\rho\vec u\otimes\vec u + p I) &= 0,\\
\partial_t E + \nabla\cdot\big(\vec u\,(E + p)\big) &= 0,
\end{aligned}
\qquad
E = \frac{p}{\gamma - 1} + \frac{1}{2}\rho\lVert\vec u\rVert^2 .
```

In a configuration: `:PDE_name => :euler, :PDE_D => D, :PDE_gamma => 1.4`. Initial states are
given as conservative variables.

## Implementation

- **Construction.** `EulerEquation(Val(D), T, gamma, rep)` sets `M = D + 2`. `rep` is the
  representation tag (default `Conservative()`).
- **Variables.** `prim2cons` maps ``(\rho, \vec u, p) \mapsto (\rho, \rho\vec u, E)``. `cons2prim`
  is its inverse and floors density and pressure at ``10^{-7}``.
- **Flux.** Defined for the `Conservative` representation.
- **Wave speed.** `max_eigenvalue` returns ``|u_d| + c`` with ``c = \sqrt{\gamma p / \rho}``, using
  a ``10^{-9}`` floor on density and pressure. For `Conservative` it expects conservative
  variables; for an `NCRepresentation` it expects primitive variables.

- **Admissibility.** [`is_admissible`](@ref) accepts finite conservative states with positive
  density and pressure. MOOD flags every candidate that fails it, whatever its criterion says
  (HyperCloudCore's [PAD](@extref HyperCloudCore Physical-admissibility-detection-(PAD))).
- **Moving grids.** `fluid_velocity` returns ``\vec u``, so the Lagrangian movers work.

## Exact solution: 1D Riemann problem

For `Riemann` initial data in 1D the closure computes the exact solution (e.g. for the Sod
shock tube):

1. **Star pressure.** ``p^*`` is found with Newton's method on the pressure function, using a
   finite-difference derivative and the tolerance ``|f(p^*)| < 10^{-9}``. This gives ``u^*``.
2. **Wave types.** Each outer wave is a shock if ``p^* > p_{L/R}``, otherwise a rarefaction.
3. **Evaluation.** At a point ``(x, t)`` the closure selects the region by the similarity
   variable ``(x - x_0)/t``: outer state, rarefaction fan, star state left or right of the
   contact. It returns conservative variables.

If the normal points in ``-x``, the problem is solved in the mirrored frame, so both
orientations of the same data give the same solution. Vacuum generation is not handled: ``p^*``
is floored at ``10^{-9}``.

## API

```@docs
EulerEquation
HyperCloudCore.is_admissible(::EulerEquation{D, M, T, <:Conservative}, ::State{M, T}) where {D, M, T}
```
