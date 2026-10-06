# Burgers' Equation

The `BurgersEquation` module implements the inviscid Burgers' equation for arbitrary spatial dimensions $D$. It describes the nonlinear propagation of waves where the wave speed is directly proportional to the state itself, naturally leading to the formation of shock waves and rarefaction fans.

## Mathematical Formulation

The scalar conservation law for the multi-dimensional inviscid Burgers' equation is given by:

$$\frac{\partial u}{\partial t} + \sum_{d=1}^D \frac{\partial}{\partial x_d} \left( \frac{1}{2} u^2 \right) = 0$$

## API Implementation Details

*   **State & Dimensionality:** The system enforces a strictly scalar state vector ($M=1$) across any user-defined spatial dimension $D$.
*   **Variable Conversions:** Because the system evaluates a single conserved scalar, `prim2cons` and `cons2prim` default to an identity mapping (`u -> u`).
*   **Flux Evaluation:** The physical flux $F_d(u)$ computes the quadratic term $0.5 u^2$ across all spatial dimensions.
*   **Eigenvalues:** The maximum eigenvalue (wave speed) along any axis simplifies to the absolute local state value, $\vert{}u\vert{}$.
*   **Velocity (Jacobian):** The advective velocity matrix used for upwind evaluations is exactly the local state $u$, returned as a strictly typed $1 \times 1$ static matrix.

## Analytic Closures

The module provides extensive support for exact analytical solutions, known as closures, which evaluate the exact mathematical state at any continuous spacetime coordinate $(x, t)$ without requiring numerical integration.

*   **1D Riemann & Smoothed Riemann (`Riemann`, `SRiemann`):** Computes exact physical shock speeds via the Rankine-Hugoniot condition and traces linear expansion fans for rarefactions.
*   **1D Box (`Box`):** Solves exact wave interactions for both "top-hat" (producing trailing shocks) and "well" (producing trailing rarefactions) initial conditions, accurately predicting shock-rarefaction intersection times.
*   **1D Sine (`Sine`):** Resolves the smooth non-linear wave steepening over time using a Newton-Raphson characteristic trace.
*   **Multi-Dimensional Gauss (`Gauss`):** Employs an iterative Newton-Raphson root-finding algorithm to trace back characteristics for Gaussian pulses in N-dimensions.
*   **2D Riemann (`Riemann`):** Extends the 1D shock and rarefaction fan logic using the directional normal vector of the interface.

## Documentation

```@docs
BurgersEquation
```