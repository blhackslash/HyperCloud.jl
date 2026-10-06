# Linear Advection System

The `LinearAdvection` module implements N-dimensional, M-variable coupled linear advection systems. It models the passive transport of state variables along predefined constant velocity fields.

## Mathematical Formulation

The coupled linear advection system is defined by a set of constant velocity matrices $A_d$ for each spatial dimension $d$:

$$\frac{\partial U}{\partial t} + \sum_{d=1}^D A_d \frac{\partial U}{\partial x_d} = 0$$

where $U$ is the state vector of length $M$, and $A_d$ is an $M \times M$ matrix defining the advection speeds and inter-variable couplings along the $d$-th axis.

## API Implementation Details

*   **Initialization & Precomputation:** The `LinearAdvection` constructor automatically extracts the spatial dimension $D$ and system size $M$ from the provided velocity matrices. It precomputes the exact spectral radius (the maximum absolute eigenvalue) for each spatial matrix using `maximum(abs.(eigvals(Matrix(vel_svec[d]))))` to ensure highly optimized runtime execution.
*   **Variable Conversions:** Linear advection does not differentiate between primitive and conservative states; `prim2cons` and `cons2prim` apply an identity mapping.
*   **Flux Evaluation:** The flux simply evaluates the matrix-vector product $A_d U$ for each spatial dimension.
*   **Eigenvalues & Velocity:** `max_eigenvalue` instantly returns the precomputed spectral radius, while `velocity` returns the stored $M \times M$ constant Jacobian $A_d$.
*   **Kinetic Relaxation:** Exposes a `kinetic_wave_speed` method that targets the specific diagonal element of the velocity matrix, supporting specialized kinetic relaxation solvers.

## Analytic Closures

The module includes a universal analytic closure capable of resolving the exact spacetime state of any coupled linear system, subject to physical validity constraints.

*   **Commutativity Verification:** For multi-dimensional domains ($D > 1$), the exact solution strictly requires that the velocity matrices commute ($A_i A_j = A_j A_i$). The closure automatically validates this via an approximate equality check (`isapprox(..., atol=1e-12)`) and throws an explicit error if the system cannot be solved via multi-dimensional ray-tracing.
*   **Eigendecomposition:** The solver computes the left and right eigenvectors ($L$ and $R$) of the summed velocity matrix, enabling a clean transformation into decoupled characteristic variables.
*   **Characteristic Tracing:** To evaluate a point at time $t$, the closure traces the characteristic waves backward to the origin $x_0 = x - \lambda t$, explicitly wrapping the coordinates if periodic boundary conditions are active. It queries the user's defined initial condition at $x_0$, projects it into characteristic space via $L$, and reconstructs the physical state via $R$.

## Documentation

```@docs
LinearAdvection