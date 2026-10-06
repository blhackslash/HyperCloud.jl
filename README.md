# HyperCloud
**HyperCloud.jl** is an advanced meshfree simulation framework designed to solve arbitrary hyperbolic conservation laws. Powered by a high-performance, physics-agnostic mathematical engine, HyperCloud provides a complete computational environment. It combines dynamic particle topologies with high-order spatial approximations—including moving-least-squares MUSCL reconstruction and Multi-Dimensional Optimal Order Detection (MOOD)—controlled entirely through a streamlined, macro-driven execution workflow.

HyperCloud integrates natively with the **PDEStudio** ecosystem to deliver a frictionless research pipeline from computation to visualization. All simulation outputs are rigorously managed and serialized via `PDEStudioCore.jl`, providing reproducible, cryptographically hashed disk caching and highly optimized statistical integration. For visual analysis, HyperCloud interfaces directly with `PDEStudio.jl`, allowing you to instantly route your meshfree Lagrangian datasets to a high-performance local renderer (`GLMakie`) or serve interactive web visualizations directly from a headless compute node to your browser (`WGLMakie` via Bonito).

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