# Kinetic Relaxation Source Terms

The kinetic relaxation framework provides the mathematical and programmatic infrastructure to approximate complex macroscopic conservation laws using simplified, linearly-advected kinetic equations. The system is driven toward the macroscopic physical state via a stiff relaxation source term.

## Kinematic-Macroscopic Mapping (`Kin2Macro`)

The `Kin2Macro{NM, NK}` structure acts as the structural bridge between the $NK$-dimensional kinetic state and the $NM$-dimensional macroscopic state. 

*   **Initialization:** It is constructed using a vector or tuple of edge indices, generating exact `ranges` for the kinetic variables that correspond to each macroscopic variable, alongside a direct `k_to_m` lookup array.
*   **State Reduction:** When applied to a kinetic vector, the functor automatically reduces it to the macroscopic state by summing the corresponding kinetic components.
*   **Index Lookup:** When applied to an integer $k$, it returns the macroscopic index $m$ associated with that kinetic dimension.

## Local Relaxation Source Term (`RelaxationSourceTerm`)

The `RelaxationSourceTerm` models the local collision-like relaxation of the kinetic system toward the macroscopic equilibrium (the Maxwellian).

### Precomputation & Setup
To maximize runtime performance, the constructor precomputes the `scaled_inv_speeds` for all $NK$ kinetic components. It queries the exact kinetic wave speeds via the `kinetic_wave_speed` API, scaling them by an interior factor (defaulting to the spatial dimension $D$) and storing their inverses. It also stores the inverse of the relaxation parameter (`inv_epsilon`) to replace costly divisions with multiplications during the time-stepping loop.

### Mathematical Formulation
The macroscopic equilibrium state (the Maxwellian) for the $k$-th kinetic variable is evaluated as:

$$M_k = c_m \left( u_m + \sum_{d=1}^D F_d^m \lambda_{d, k}^{-1} \right)$$

where $c_m$ is the base coefficient for the associated macroscopic variable $m$, $u_m$ is the macroscopic state, $F_d^m$ is the physical macroscopic flux, and $\lambda_{d, k}^{-1}$ is the precomputed scaled inverse kinetic speed.

The continuous source term evaluated by `evaluate_source` is then simply the linear relaxation:

$$S_k = \frac{M_k - U_k}{\epsilon}$$

where $U_k$ is the current kinetic state.

### Implicit Time Integration
Because the relaxation parameter $\epsilon$ can be extremely small (approaching the zero-relaxation limit), explicit evaluation of this source term would impose severe stability restrictions on the timestep. To bypass this, the framework implements an exact `implicit_solve` function.

The implicit update analytically resolves the stiff differential equation within the Runge-Kutta stage, updating the kinetic state algebraically:

$$U_k^{n+1} = \frac{Y_{in, k} + \frac{\Delta t}{\epsilon} M_k}{1 + \frac{\Delta t}{\epsilon}}$$

where $Y_{in}$ is the incoming kinetic state for the current stage, and $\Delta t$ is the scaled timestep coefficient.

## Kinetic Initialization

To seamlessly integrate with the broader HyperCloud macroscopic workflow, the framework overloads `set_initial_conditions!` for kinetic source terms. 

Instead of requiring the user to define complex $NK$-dimensional initial conditions, the user provides a standard macroscopic `InitialCondition`. The function iterates over the `ParticleGrid`, evaluates the macroscopic state and physical fluxes at each coordinate, and directly initializes the particles into the kinetic Maxwellian equilibrium state.

## Documentation

```@docs
Kin2Macro
RelaxationSourceTerm
```