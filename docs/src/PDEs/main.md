# PDE Architecture & Exact Solutions

Every equation is a struct subtyping `HyperbolicPDE{D, M, T}`, where `D` is the spatial
dimension, `M` the number of variables and `T` the floating-point type. The numerical schemes
only access an equation through a small set of functions, so new physics plugs in without
touching the solver.

## Equation interface

| Function | Required | Returns |
|:---------|:---------|:--------|
| `flux(eq, U)` | yes | `Flux{D, M, T}`: the physical flux ``F_d(U)`` for every direction `d`. |
| `max_eigenvalue(eq, U, d)` | yes | Largest absolute wave speed in direction `d` (used for the CFL step and Rusanov flux). |
| `prim2cons(eq, V)`, `cons2prim(eq, U)` | yes | Conversion between primitive and conservative variables (identity for scalar laws). |
| `velocity(eq, U, d)` | upwind schemes | Flux Jacobian ``\partial F_d / \partial U`` as an `M × M` static matrix. |
| `kinetic_wave_speed(eq, d, k)` | kinetic systems | Speed of kinetic variable `k` in direction `d`. |
| `analytic_closure(eq, ic, geom)` | optional | Exact solution, see below. |

A builder method makes the equation available from configurations:

```julia
function HyperCloud.build_equation(::Val{:my_pde}, pde_conf::Dict, context::Dict)
    T = context[:T]::DataType
    return MyPDE{pde_conf[:D], T}()
end
```

## Representations

Some equations exist in several formulations. They carry a representation tag, a subtype of
`EquationRepresentation`. `Conservative` is the conservative form; `NCRepresentation` covers
non-conservative (e.g. primitive) forms, whose integration paths are configured through
`parse_representation`. Methods such as `flux` and `max_eigenvalue` dispatch on the tag.

## Exact solutions

[`analytical_solution`](@ref) builds the equation, domain and initial condition from a
configuration and calls [`analytic_closure`](@ref). The result is a closure

```julia
exact = analytical_solution(params)
exact(SVector(x..., t))   # -> State{M}
```

that evaluates the exact solution at any spacetime point. It is used as the reference for
error statistics and plots. The reference is computed in `params[:real_type]` (default
`Float64`), independently of the simulation precision.

Exact solutions currently available:

| Equation | Initial condition | Dimension |
|:---------|:------------------|:----------|
| `LinearAdvection` | any | any `D` (commuting velocity matrices) |
| `BurgersEquation` | `Riemann`, `SRiemann` | any `D` |
| `BurgersEquation` | `Gauss` | any `D`, before the breaking time |
| `BurgersEquation` | `Box` | 1D |
| `EulerEquation` | `Riemann` | 1D |

Other combinations fall back to a warning. Burgers returns the initial condition (exact only
at `t = 0`); every other equation returns `nothing`.

To add an exact solution, add a method for your combination:

```julia
function HyperCloud.analytic_closure(eq::MyPDE{D, T}, ic::Gauss, geom::GeometricDomain) where {D, T}
    # precompute everything that does not depend on (x, t) here
    return function (st::SVector)
        t = T(st[end])
        x = SVector{D, T}(ntuple(d -> T(st[d]), Val(D)))
        # ...
        return State{1, T}(u)
    end
end
```

`geom` gives access to the domain bounds (`geom.mins`, `geom.maxs`) and periodicity
(`geom.is_periodic`).

## API

```@docs
analytical_solution
analytic_closure
EquationRepresentation
```
