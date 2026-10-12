# Linear Advection

`LinearAdvection` implements coupled linear transport with constant matrices:

```math
\frac{\partial U}{\partial t} + \sum_{d=1}^D A_d \frac{\partial U}{\partial x_d} = 0,
```

where ``U`` has ``M`` components and every ``A_d`` is an ``M \times M`` matrix.

## Specifying the velocities

`velocities` has one entry per dimension, each describing ``A_d``. The dimension `D` is the
number of entries.

```julia
LinearAdvection(((1.0,),), Float64)               # 1D scalar: u_t + u_x = 0
LinearAdvection(((1.0,), (0.5,)), Float64)        # 2D scalar, velocity (1, 0.5)
LinearAdvection(((1.0, -1.0),), Float64)          # 1D, two variables: diag(1, -1)
```

A tuple of ``M`` numbers is read as the diagonal matrix ``\operatorname{diag}(a_1, \dots, a_M)``.
This is the format produced by [`generate_relax_params`](@ref) for kinetic systems. Other
formats are converted by `param2vel`.

In a configuration: `:PDE_name => :linear, :PDE_velocities => ((1.0,), (1.0,))`.

## Implementation

- **Wave speeds.** The constructor precomputes the spectral radius of every ``A_d``;
  `max_eigenvalue` returns it without any work at run time.
- **Variables.** `prim2cons` and `cons2prim` are the identity.
- **Flux.** ``F_d(U) = A_d U``; `velocity(eq, U, d)` returns ``A_d``.
- **Kinetic use.** `kinetic_wave_speed(eq, d, k)` returns the diagonal entry
  ``(A_d)_{kk}``. It is meaningful for diagonal systems, which is how kinetic relaxation
  systems are built.

## Exact solution

The exact solution works with any initial condition and any dimension:

1. **Commutation check.** For ``D > 1`` the matrices must commute
   (``A_i A_j = A_j A_i``). Otherwise the closure throws an error.
2. **Characteristic variables.** It diagonalizes ``\sum_d A_d = R \Lambda R^{-1}``. The speed of
   characteristic ``m`` in direction ``d`` is ``\lambda_{m,d} = (R^{-1} A_d R)_{mm}``.
3. **Tracing.** Each characteristic variable is the initial condition at
   ``x - \lambda_m t``, projected with ``R^{-1}``. Periodic directions are wrapped into the domain.
   The physical state is reassembled with ``R``.

!!! note
    Step 2 assumes that the eigenvectors of ``\sum_d A_d`` also diagonalize every ``A_d``. This
    always holds for diagonal systems. It can fail for non-diagonal matrices when
    ``\sum_d A_d`` has repeated eigenvalues.

## API

```@docs
LinearAdvection
```
