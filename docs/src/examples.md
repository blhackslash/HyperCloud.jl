# Included Examples & Quick Start

HyperCloud ships with a suite of standard benchmark configurations demonstrating the framework's capabilities across different physical systems, grid topologies, and numerical schemes. 

## Quick Start Tutorial

To run any of the built-in examples, you only need to load the framework alongside a Makie backend and the `PDEStudio` visualization suite. 

To execute the simulation and instantly visualize the results in an interactive window, use the `@plot` macro with the uncalled example function name:

```julia
using HyperCloud
using PDEStudio
using GLMakie # Or WGLMakie for browser-based rendering

# Launch the interactive plotter for the 1D Sod Shock Tube
@plot euler_1d_sod
```

If you want to run the simulation headlessly without visualization (e.g., for benchmarking or cluster execution), use the `@sim` macro instead:

```julia
using HyperCloud

# Execute the simulation headlessly
@sim euler_1d_sod
```

## Available Benchmarks

### 1D Sod Shock Tube (`euler_1d_sod`)
A classic 1D Riemann problem for the Euler equations. 
*   **Domain:** Simulates an interface separating a high-pressure left state ($\rho=1.0, u=0.0, p=1.0$) and a low-pressure right state ($\rho=0.125, u=0.0, p=0.1$).
*   **Grid:** Resolves the domain using 330 particles bounded by a Fixed Dirichlet condition on the left and an Outflow condition on the right.
*   **Schemes:** Compares the exact analytical solution against a 1st-order baseline (`RK2Upwind1`), a TVD limiter (`RK2MUSCL2_VK`), and high-order MOOD strategies (e.g., `RK4MUSCL4_U2_EPD1`).

### 1D Kinetic Burgers (`burgers_kinetic`)
Demonstrates the integration of the kinetic relaxation solver applied to the 1D Burgers' equation.
*   **Domain:** Advects a Gaussian pulse across a periodic domain using 300 particles.
*   **Kinetic Setup:** Generates symmetric kinetic velocities with a magnitude of 2.0 and enforces a stiff kinetic relaxation parameter ($\epsilon = 10^{-6}$).
*   **Schemes:** Directly compares a standard macroscopic explicit scheme (`RK2MUSCL2`) against the stiff kinetic IMEX formulation (`ARS222MUSCL2`).
*   **Parameter Sweeps:** Uses a varied parameter dictionary to sweep the execution across different grid seeds (`10, 100, 1000, 10000`).

### 2D Linear Advection (`advection_2d`)
A baseline multi-dimensional transport test.
*   **Domain:** Advects a Gaussian pulse across a 2D periodic rectangular domain using 100x100 particles.
*   **Schemes:** Evaluates the `RK2MUSCL2` spatial reconstruction against the exact characteristic ray-tracing analytical solution.
*   **Parameter Sweeps:** Evaluates performance and stability across varied grid randomization seeds (`10, 100, 1000, 10000`).
