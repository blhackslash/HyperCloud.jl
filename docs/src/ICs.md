# Initial Conditions

The initial condition architecture provides a library of composable functors that evaluate the starting state of the physical system at any given spatial coordinate.

## Core Architecture

Initial conditions in HyperCloud are implemented as callable structs (functors) derived from the abstract type `InitialCondition`. Each functor accepts an N-dimensional position vector (`Space{D}`) and returns a strictly typed physical state vector (`State{M}`). 

The application of these conditions is handled universally by `set_initial_conditions!`, which iterates over the particle grid and evaluates the designated functor for each particle's spatial position. The instantiation of these structs is managed by the `build_ic` factory function, which uses dynamic `Val` dispatch based on the `:name` key in the configuration dictionary to construct the appropriate geometry.

## Available Geometries

*   **Box (`Box`):** Defines a discontinuous rectangular sub-region. Returns `u_box` if the queried position falls strictly within the defined `mins` and `maxs` boundaries, and `u_bg` otherwise.
*   **Gaussian Pulse (`Gauss`):** Evaluates a smooth, N-dimensional Gaussian bell curve. The state is computed as $a \exp(-\frac{\Vert{}pos - b\Vert{}^2}{width^2})$, where $a$ is the amplitude state, $b$ is the spatial center, and $width$ controls the spread.
*   **Standard Riemann (`Riemann`):** Represents a classic shock tube or half-space problem. It splits the domain using a hyperplane defined by an origin point `p0` and a normal vector `n`. Returns the left state `uL` if the dot product of the relative position and the normal is negative, and the right state `uR` otherwise.
*   **Smoothed Riemann (`SRiemann`):** A continuous approximation of the standard Riemann problem. It utilizes an arctangent function to smoothly transition between `uL` and `uR` over a specified `width` parameter, preventing infinite gradients at the interface.
*   **Quadrant Riemann (`QRiemann`):** Partitions the space into $2^D$ distinct orthogonal regions intersecting at a central point `p0`. It uses a binary encoding strategy based on the spatial axes to map the evaluated coordinate to one of the provided $N\_states$.
*   **Sinusoidal Wave (`Sine`):** Evaluates a smooth, multidimensional periodic wave. The physical state is calculated as $a \sin(2\pi \sum \frac{pos}{period}) + c\_offset$.

## API Documentation

```@docs
InitialCondition
set_initial_conditions!
Box
Gauss
QRiemann
Riemann
Sine
SRiemann
```