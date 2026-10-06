# PDE Architecture & Analytical Solutions

The framework utilizes a highly modular architecture to parse physical configurations, construct equation systems, and generate exact mathematical closures. This decoupled design ensures that numerical solvers operate independently of the underlying physics.

## Equation Representations & Paths

The framework uses a type-hierarchy to handle different formulations of the governing equations.

*   **Representations:** Equations are categorized via abstract types like `EquationRepresentation` and `NCRepresentation`. Concrete types like `Conservative` dictate how the flux and eigenvalues are mathematically evaluated by the solver.
*   **Non-Conservative Paths:** For primitive or Lagrangian formulations, the configuration dictionary parses specific integration path topologies, resolving symbols to concrete types such as `MappedPath`, `LinePath`, or `NaiveAveragePath`.

## The Equation Builder

All PDE structs are instantiated through a unified factory pattern.

*   **`build_equation`:** This function requires the configuration dictionary to include a strictly typed `:name` symbol (e.g., `:linear` or `:burgers`). It dynamically dispatches to the correct specific constructor using `Val(eq_name)`, passing along the isolated PDE parameters and the execution context.

## Analytical Solution Generator

To support robust verification and error analysis, the framework provides a unified pipeline for generating exact mathematical solutions.

*   **`analytical_solution`:** This orchestrator function takes the global shared parameter dictionary and isolates the nested `PDE`, `Grid`, and `IC` configurations. 
*   **Assembly:** It sequentially builds the PDE struct, instantiates the continuous `GeometricDomain` (which encapsulates boundary and periodicity logic), and constructs the `InitialCondition`.
*   **Closure Generation:** Finally, it dispatches these three core components to `analytic_closure`, which returns a fast, standalone closure: `exact_u(st)`. This resulting function accepts a unified spacetime tensor and evaluates the exact analytical state at any continuous coordinate $(x, t)$.

## Documentation

```@docs
analytical_solution
analytic_closure
EquationRepresentation
```