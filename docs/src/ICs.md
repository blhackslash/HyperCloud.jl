# Initial Conditions

An initial condition is a callable struct that subtypes `InitialCondition`. It takes a position
`pos::Space{D}` and returns the state `State{M}` at that point.

[`set_initial_conditions!`](@ref) evaluates it at every particle. The kinetic pipeline uses a
second method that maps the macroscopic state to the kinetic equilibrium (see
[Kinetic Relaxation](SourceTerms/kinetic.md)).

In a configuration the initial condition is selected with `:IC_name`; [`build_ic`](@ref)
converts all `IC_` parameters to the simulation type `T`. The keys of each type are listed in
the [Configuration Reference](configuration.md).

## Available initial conditions

| Type | `:IC_name` | State at `pos` |
|:-----|:-----------|:---------------|
| [`Gauss`](@ref) | `:gauss` | ``a \exp\left(-\lVert pos - b \rVert^2 / \mathrm{width}^2\right)`` |
| [`Box`](@ref) | `:box` | `u_box` if `mins .<= pos .<= maxs` in every direction, else `u_bg` |
| [`Sine`](@ref) | `:sine` | ``a \sin\left(2\pi \sum_d pos_d / \mathrm{period}_d\right) + c_\mathrm{offset}`` |
| [`Riemann`](@ref) | `:riemann` | `uL` if ``(pos - p_0) \cdot n < 0``, else `uR` |
| [`SRiemann`](@ref) | `:s_riemann` | ``\frac{u_L + u_R}{2} - \frac{u_L - u_R}{\pi} \arctan\left(\frac{(pos - p_0)\cdot n}{\mathrm{width}}\right)`` |
| [`QRiemann`](@ref) | `:q_riemann` | `u_states[1 + Σ_d 2^(d-1) (pos_d ≥ p0_d)]` |

Notes:

- **Riemann normals.** `build_ic` normalizes `n` (default ``+x``). `uL` lies on the side the
  normal points away from.
- **Quadrant ordering.** `QRiemann` needs `2^D` states in binary order. In 2D:
  `(bottom-left, bottom-right, top-left, top-right)`.
- **Precision.** `Gauss` and `SRiemann` store their parameters in the precision of the
  amplitude/state. A Float32 configuration therefore evaluates entirely in Float32.

## Adding an initial condition

Define a struct, make it callable, and add a builder method:

```julia
using HyperCloud, StaticArrays

struct Ring{D, M, T} <: InitialCondition
    u_in::State{M, T}
    u_out::State{M, T}
    center::Space{D, T}
    radius::T
end

(ic::Ring)(pos::Space{D}) where {D} = sum(abs2, pos - ic.center) < ic.radius^2 ? ic.u_in : ic.u_out

function HyperCloud.build_ic(::Val{:ring}, conf::Dict, ctx::Dict)
    T, D, M = ctx[:T]::DataType, ctx[:D]::Int, ctx[:M]::Int
    return Ring(State{M, T}(Tuple(T.(conf[:u_in]))), State{M, T}(Tuple(T.(conf[:u_out]))),
                Space{D, T}(Tuple(T.(conf[:center]))), T(conf[:radius]))
end
```

`:IC_name => :ring` now selects it. To use it as a reference solution, add an
[`analytic_closure`](@ref) method for the equations where it has one (see
[PDE Architecture](PDEs/main.md)).

## API

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
