# Examples

HyperCloud ships with benchmark configurations covering the built-in equations and the simulation
pipelines.

## Running an example

Pass the example function, uncalled, to a macro:

```julia
using HyperCloud
using PDEStudio
using GLMakie          # or WGLMakie for browser-based rendering

@plot euler_1d_sod     # run (or load from cache) and open the viewer
```

Without a viewer, for example on a cluster:

```julia
using HyperCloud

@sim euler_1d_sod
```

## Available benchmarks

### 1D Sod shock tube (`euler_1d_sod`)

The classic Riemann problem for the Euler equations with ``\gamma = 1.4``.

- **Data:** ``(\rho, u, p) = (1, 0, 1)`` for ``x < 0`` and ``(0.125, 0, 0.1)`` for ``x > 0``,
  given to the IC as conservative states. `SOD_LEFT` and `SOD_RIGHT` hold them.
- **Grid:** 330 jittered particles on ``[-0.5, 0.5]``, fixed Dirichlet on the left and outflow
  on the right, `tmax = 0.2`.
- **Methods:** the exact solution, first-order `RK2Upwind1`, `RK2MUSCL2_VK` (Venkatakrishnan
  limiter), and MUSCL with MOOD (`RK2MUSCL2_U2_EPD1`, `RK4MUSCL4_U2_EPD1`, plus the inactive
  `RK2MUSCL2_U2_EPD2` and `RK2MUSCL2_U1_EPD1`).

### 1D kinetic Burgers (`burgers_kinetic`)

Compares a direct scheme with the kinetic relaxation formulation for Burgers' equation.

- **Data:** a Gaussian pulse (`a = 1`, `width = 1`) on the periodic domain ``[-5, 5]`` with
  300 particles, up to `tmax = 1.0`. The pulse breaks at ``t^* \approx 1.17``, so the run stays
  in the smooth regime where the exact solution is available.
- **Kinetic setup:** speeds ``\pm 2`` from `generate_relax_params(2.0, 1, 1)` and
  ``\epsilon = 10^{-6}``.
- **Methods:** `RK2MUSCL2` (direct) and `ARS222MUSCL2` (kinetic, IMEX).
- **Sweep:** grid seeds `10, 100, 1000, 10000`.

### 1D BGK shock tube (`bgk_1d_sod`)

The Sod problem for the discrete-velocity [BGK model](PDEs/bgk.md) with ``\gamma = 1.4``.

- **Model:** 24 uniform velocity nodes on ``[-6, 6]`` and the Chu-reduced pair ``(f, h)``
  (48 equations); the IC gives the Sod states as Euler variables.
- **Grid:** as `euler_1d_sod` (330 jittered particles, Dirichlet left, outflow right, `tmax = 0.2`).
- **Methods:** `ARS222MUSCL2_VK` and `ARS222MUSCL2_U2_EPD1` (implicit collision), the inactive explicit
  `RK2MUSCL2_VK` (only stable for the larger ``\tau``), against the exact Euler solution
  (`reference = :euler_limit_solution`).
- **Sweep:** relaxation times ``\tau = 10^{-5}, 10^{-3}, 10^{-2}``: the hydrodynamic limit, then
  increasingly viscous and heat-conducting profiles.

### 2D linear advection (`advection_2d`)

Transport of a Gaussian pulse along the diagonal.

- **Data:** velocity ``(1, 1)`` on the periodic square ``[-5, 5]^2`` with ``50 \times 50``
  jittered particles, up to `tmax = 10` (one full period in each direction).
- **Methods:** `RK2MUSCL2` against the exact solution from characteristic tracing.
- **Sweep:** grid seeds `10, 100, 1000, 10000`.

## Writing your own

Copy one of the example functions, change the parameters, and run it with `@sim` or `@plot`.
The [Configuration Reference](configuration.md) lists every key. Saving the function as
`my_case.jl` in a folder listed in `CONFIG_SEARCH_PATHS` lets you call `@sim my_case` from any
session without including the file first.
