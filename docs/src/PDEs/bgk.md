# BGK Equation

[`BGKEquation`](@ref) is a discrete-velocity model of the Bhatnagar–Gross–Krook equation

```math
\partial_t f + \vec v\cdot\nabla_x f = \frac{1}{\tau}\big(\mathcal M[f] - f\big),
```

for the distribution ``f(\vec x, \vec v, t)`` of a gas in ``D`` space dimensions. Its moments
``U = (\rho, \rho\vec u, E)`` are the conservative variables of the Euler equations, which the BGK
solution approaches for ``\tau \to 0``; for ``\tau \to \infty`` every velocity moves freely.

## Velocity space as the state vector

The velocity space is discretised by a tensor grid of `Nv` nodes per direction
([`velocity_grid`](@ref)), and the particle state is the vector of distribution values at the
nodes: the number of equations is the number of discrete velocities. Everything HyperCloudCore
offers (MUSCL, limiters, MOOD, IMEX, moving grids, the device backends) therefore runs on the BGK
model unchanged. Moments are quadrature sums,
``U = \sum_k w_k\,(1, \vec v_k, \tfrac12|\vec v_k|^2)\,f_k``.

Since the state is an `SVector`, its length drives compile time and per-particle cost. The model
keeps it small in three ways:

- **Only `D` velocity directions.** Degrees of freedom beyond the resolved velocity directions (the
  remaining directions of a 3-D gas and internal ones, ``K = 2/(\gamma - 1) - D`` in total) are
  carried by a second distribution ``h = \int \tfrac12|\xi|^2 f\,d\xi`` on the same nodes (a Chu
  reduction). The state is ``f`` for ``\gamma = (D + 2)/D`` (``M = N_v^D``), otherwise ``(f, h)``
  (``M = 2 N_v^D``). A 1-D shock tube with ``\gamma = 1.4`` needs a 1-D velocity grid, not ``N_v^3``
  velocities.
- **Gauss–Hermite nodes** (`:velocity_grid => :gauss_hermite`) integrate polynomials times the
  Gaussian of their scale exactly, so near-equilibrium flows need far fewer nodes than the uniform
  midpoint grid (`:uniform`, the robust choice for strongly non-equilibrium flows).
- **No per-velocity arrays in kernels.** The equation is a small bits type (1-D nodes and weights,
  grid centre, ``K``, ``\gamma``, relaxation model); the node of a component is resolved from its index
  at compile time, so the flux and the dissipation unroll into straight-line code.

## Collision term

The equilibrium is the conservative discrete Maxwellian of [`equilibrium`](@ref)
(``f_k = \exp(\alpha\cdot\psi_k)`` with ``\psi = (1, \vec v, \tfrac12|\vec v|^2)``, Mieussens'
discrete entropic equilibrium, ``h_k = \tfrac K2\,\hat\theta f_k``): ``\alpha`` is found by a damped
Newton iteration on the convex dual problem, started from the continuous Maxwellian, so the discrete
moments of the equilibrium equal ``U`` to round-off and the collision term conserves mass, momentum
and energy exactly. A few microseconds per particle; the iteration runs inside the kernels.

The relaxation time comes from a [`RelaxationTimeModel`](@ref): constant ([`ConstantRelaxationTime`](@ref),
`Inf` for free transport) or ``\tau = \mu(\theta)/p`` with a power-law viscosity
([`ViscosityRelaxationTime`](@ref)).

With an IMEX tableau the collision is the implicit source [`BGKRelaxation`](@ref). Because it
conserves the moments, the implicit stage equation has the closed-form solution
``f = (f^* + \lambda\,\mathcal M[f^*])/(1 + \lambda)``, ``\lambda = \Delta t\,a_{ii}/\tau``, for any
stiffness, which makes the scheme asymptotic preserving: at ``\tau = 10^{-5}`` the BGK Sod solution
is as close to the exact Euler solution as the Euler scheme itself. With an explicit tableau the
collision becomes [`ExplicitBGKRelaxation`](@ref) (stable for ``\Delta t \lesssim \tau`` only).

## Transport

The flux is ``v_{k,d} f_k``. The transport is linear and diagonal, so the Rusanov flux of every
velocity is its exact upwind flux: `RusanovFlux`, `UpwindFlux` and `AlignedRusanovFlux` use the
dissipation ``|v_{k,d}|`` per component instead of the largest node speed for all components (which
smears the slow velocities; on a smooth collisionless wave it costs a factor 2 to 10 in the error).
`max_eigenvalue` is the largest node speed, which sets the CFL time step.

The PAD of MOOD (`is_admissible`) asks for positive density and temperature; individual ``f_k`` may
dip below zero.

## Configuration

```julia
:PDE_name => :bgk, :PDE_D => 1, :PDE_gamma => 1.4,
:PDE_Nv => 24, :PDE_vmax => 6.0, :PDE_velocity_grid => :uniform,
:PDE_tau => 1e-3,                       # or :PDE_viscosity => μ_ref, :PDE_omega => ω
:sim_func_name => :run_bgk_simulation,  # direct and kinetic runners forward BGK runs here
```

Initial conditions give the macroscopic state ``(\rho, \rho\vec u, E)``, exactly as for an Euler
run with the same ``\gamma``; every particle starts in the discrete Maxwellian of that state. The
stored snapshots are the moments, so BGK and Euler runs can be compared directly
(`:save_relax => true` keeps the distributions). Choose `vmax` so that the Maxwellians of the flow
decay inside the grid (``|u| + 5\sqrt\theta`` or more) and enough nodes to resolve them
(``\Delta v \lesssim \sqrt\theta``).

References: the collisionless model has the exact solution of [`analytic_closure`](@ref)
(each velocity transported, initial data in equilibrium); for small ``\tau`` use
[`euler_limit_solution`](@ref) as the `reference` of a `SimulationConfig` (see the example
`bgk_1d_sod`).

## Cost of the velocity space

Measured on 8 threads (ARS222, MUSCL3 with the Venkatakrishnan limiter, regular lattice, `Float64`),
one time step and the first (compiling) step:

| D | M | particles | step / particle | first step | per-particle buffers | neighbour-slot buffers |
|:--|--:|--:|--:|--:|--:|--:|
| 1 | 16 | 4000 | 0.32 µs | 3 s | 1.9 KiB | 10 × 0.27 KiB |
| 1 | 64 | 4000 | 1.6 µs | 5 s | 6.8 KiB | 10 × 1.0 KiB |
| 1 | 128 | 4000 | 3.2 µs | 11 s | 13 KiB | 10 × 2.0 KiB |
| 1 | 256 | 4000 | 7.1 µs | 30 s | 26 KiB | 10 × 4.0 KiB |
| 2 | 16 | 64² | 2.4 µs | 5 s | 3.1 KiB | 36 × 0.41 KiB |
| 2 | 64 | 64² | 5.8 µs | 15 s | 9.1 KiB | 36 × 1.5 KiB |
| 2 | 144 | 64² | 16 µs | 40 s | 19 KiB | 36 × 3.4 KiB |
| 2 | 256 | 64² | 30 µs | 95 s | 33 KiB | 36 × 6.0 KiB |

The cost per particle grows linearly with ``M`` up to a few hundred equations. In 2-D the per-slot
buffers of the scheme (the flux differences of every neighbour, ``(D + 1) M`` values per slot)
dominate the memory, about 6× the per-particle buffers, and the compile time of the first step grows
faster than ``M``; ``16 \times 16`` velocities still work (95 s to compile, 30 µs per particle and
step). About a third of a step is the collision (two equilibria per stage), a third the gradient
pre-pass with the limiter, and the rest the flux pass.

## API

```@docs
BGKEquation
HyperCloudCore.is_admissible(::BGKEquation{D, M, T}, ::State{M, T}) where {D, M, T}
velocity_grid
num_velocities
velocity_node
moments
equilibrium
euler_limit
RelaxationTimeModel
ConstantRelaxationTime
ViscosityRelaxationTime
relaxation_time
BGKRelaxation
ExplicitBGKRelaxation
run_bgk_simulation
euler_limit_solution
bgk_1d_sod
```
