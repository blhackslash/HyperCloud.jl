# HyperCloud.jl

Physics and workflow layer on top of HyperCloudCore.jl: equations (linear advection,
Burgers, Euler), initial conditions, kinetic relaxation source terms, config-driven
builders, the simulation runners and the `@sim`/`@plot`/`@backup` macros. Reexports
HyperCloudCore and PDEStudioCore; plotting goes through PDEStudio (optional).

## Commands

```bash
julia --project -e 'using Pkg; Pkg.test()'                         # full suite (GUI tests need the X display, which the container has)
julia --project=test -e 'using Pkg; Pkg.instantiate()'              # once
julia -t 8 --project=test test/runtests.jl                          # same suite, faster to iterate
HYPERCLOUD_GUI_TESTS=false julia -t 8 --project=test test/runtests.jl   # skip PDEStudio/GLMakie tests
HYPERCLOUD_BACKEND=ka_cpu julia -t 8 --project=test test/runtests.jl    # KernelAbstractions CPU backend (checks kernel rules)
julia --project=docs docs/make.jl                                   # docs (instantiate docs env first)
julia -t auto --project=test/gpu test/gpu/runtests.jl                # device checks, AMD GPU only (CASE=<regex>)
julia -t auto --project=test/gpu benchmark/runners.jl                # CPU vs GPU runner benchmark
```

`HYPERCLOUD_BACKEND` accepts `polyester` (default), `threads`, `ka_cpu`. Run the suite
with `ka_cpu` too when a change touches anything that ends up in a kernel.

## Layout

| Path | Contents |
|:-----|:---------|
| `src/PDEs/` | Equations and representations (`Conservative`, primitive/Lagrangian paths), `analytical_solution` |
| `src/ICs/` | Initial condition functors and `build_ic` |
| `src/SourceTerms/` | Kinetic relaxation source terms |
| `src/Builder/` | `build_equation`, kinetic system builder, reference solutions |
| `src/SimulationFunctions/` | `run_direct_simulation`, `run_kinetic_simulation` |
| `src/device.jl` | Device runs: `storage_for`, `with_backend`, host/device copies, memory limit |
| `src/inbuilt_methods.jl` | `create_method`, master method dictionary |
| `src/Examples/` | Exported example functions (`euler_1d_sod`, …) used by `@sim`/`@plot` |
| `src/macros.jl` | `@sim`, `@plot`, `@backup`, Bonito/config path settings |

## Configuration model

- A run is one flat `Dict{Symbol, Any}`. Keys are namespaced by prefix (`:PDE_name`,
  `:Grid_Ns`, `:Weight_range`, `:Scheme_order`, `:MOOD_criterion`, …);
  `extract_namespace` strips the prefix before a builder sees it. Unprefixed keys are global
  (`:tmax`, `:snapshots`, `:T`, `:sim_func_name`).
- Misspelled keys are silently ignored (e.g. `:Grid_bc_map` instead of `:Grid_bc`). When
  adding or renaming a key, update `docs/src/configuration.md` and add a test that the
  key has an effect.
- Extend by dispatch on `Val(name)`: e.g. a new IC is a functor plus
  `build_ic(::Val{:my_ic}, conf::Dict, ctx::Dict)`; equations the same way via
  `build_equation`. New equations implement the HyperCloudCore API (`flux`,
  `max_eigenvalue`, `velocity`, `is_admissible`, `prim2cons`/`cons2prim`,
  `fluid_velocity` for Lagrangian movers) and should provide `analytical_solution` where
  one exists.
- Parameter dicts are hashed by PDEStudioCore to name cache files: renaming keys or
  changing defaults that end up in the dict changes the hashes of the user's existing
  results. Call that out in the PR.

## Conventions

- Requires Julia 1.12 (`[compat]`).
- Numerics that would work for any equation belong in HyperCloudCore, not here.
- Precision follows `:T`: ICs, equations and builders must keep the requested float type
  (tested for Float32 and Float64).
- Runs must work on the GPU path: setup stays on the host, the grid/scheme/stepper move to
  the device via `context[:Storage]`; no scalar indexing of device arrays in the runners.
- The runners catch exceptions and return `nothing`; tests and scripts must check
  `sim !== nothing` and should look at the log, not only at the return value.
- Every export needs a docstring that appears in the docs (`checkdocs = :exports`).
  Docs cross-link to HyperCloudCore; for a local build set `HYPERCLOUDCORE_INVENTORY` to
  `../HyperCloudCore.jl/docs/build/objects.inv` (github.io may be unreachable here).
- New examples go into `src/Examples/` (exported function + `_main.jl` include) and
  `docs/src/examples.md`; keep their default run time short.
- Commit subjects follow the repo style: `Area: what changed`.
