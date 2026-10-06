# Euler Equations

The `EulerEquation` module implements the compressible Euler equations for ideal gases. It models the dynamics of inviscid, adiabatic fluids by conserving mass, momentum, and energy across arbitrary spatial dimensions.

## Mathematical Formulation

For an ideal gas with a specific heat ratio $\gamma$, the conservative state vector is $U = [\rho, \rho \vec{u}, E]^T$. The multi-dimensional system is defined as:

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0$$
$$\frac{\partial (\rho \vec{u})}{\partial t} + \nabla \cdot (\rho \vec{u} \otimes \vec{u} + p I) = 0$$
$$\frac{\partial E}{\partial t} + \nabla \cdot (\vec{u} (E + p)) = 0$$

where the pressure $p$ is governed by the ideal gas equation of state: $E = \frac{p}{\gamma - 1} + \frac{1}{2} \rho \Vert{}\vec{u}\Vert{}^2$.

## API Implementation Details

*   **Dimensionality:** The system adapts to any spatial dimension $D$, automatically setting the state vector length to $M = D + 2$. It stores the specific heat ratio $\gamma$ alongside an equation representation tag.
*   **Variable Conversions (`prim2cons`, `cons2prim`):** Provides strict conversions between the primitive state $(\rho, \vec{u}, p)$ and the conservative state $(\rho, \rho\vec{u}, E)$. To ensure numerical robustness during aggressive transients, `cons2prim` enforces a strict floor value of $10^{-7}$ on both density and pressure.
*   **Flux Evaluation:** The `flux` function strictly evaluates the conservative mass, momentum, and energy fluxes for each spatial dimension $d$.
*   **Eigenvalues:** The `max_eigenvalue` function computes the maximum acoustic wave speed $\vert{}u_d\vert{} + c$, where $c = \sqrt{\gamma p / \rho}$. It dispatches dynamically based on whether the system is tagged as `Conservative` or `NCRepresentation` to extract the correct velocity components, utilizing a $10^{-9}$ floor for density and pressure evaluations to prevent domain errors.

## Analytic Closures

The module features an exact analytical solver for the 1D Euler Riemann problem (e.g., the Sod Shock Tube). 

*   **Iterative Solver:** The closure employs a Newton-Raphson root-finding algorithm (capped at 100 iterations) to iteratively converge on the exact pressure within the intermediate star region ($p^*$).
*   **Wave Categorization:** Based on the computed $p^*$, the solver classifies the left and right traveling waves as either discrete shocks or continuous rarefaction fans.
*   **Spacetime Evaluation:** It analytically evaluates the density, velocity, and pressure profiles across the entire domain, explicitly calculating the states within the expansion fan tails/heads and across the contact discontinuity.

## Documentation

```@docs
EulerEquation
```