# HyperCloud.jl Documentation

**HyperCloud.jl** is an advanced meshfree simulation framework designed to solve arbitrary hyperbolic conservation laws. Powered by a high-performance, physics-agnostic mathematical engine, HyperCloud provides a complete computational environment. It combines dynamic particle topologies with high-order spatial approximations—including moving-least-squares MUSCL reconstruction and Multi-Dimensional Optimal Order Detection (MOOD)—controlled entirely through a streamlined, macro-driven execution workflow.

HyperCloud integrates natively with the **PDEStudio** ecosystem to deliver a frictionless research pipeline from computation to visualization. All simulation outputs are rigorously managed and serialized via `PDEStudioCore.jl`, providing reproducible, cryptographically hashed disk caching and highly optimized statistical integration. For visual analysis, HyperCloud interfaces directly with `PDEStudio.jl`, allowing you to instantly route your meshfree Lagrangian datasets to a high-performance local renderer (`GLMakie`) or serve interactive web visualizations directly from a headless compute node to your browser (`WGLMakie` via Bonito).

---

## Project Architecture

The framework is highly modular, separating the physical modeling from the numerical execution. 

*   **`PDEs/` & `SourceTerms/`**: House the concrete physical definitions. Users define custom mathematical equations (e.g., N-dimensional Euler, Linear Advection) and volumetric source terms here, which natively plug into the solver backend.
*   **`ICs/`**: Provides a library of composable initial condition functors, such as N-dimensional Riemann problems, Gaussian pulses, Box distributions, and sinusoidal waves.
*   **`Builder/` & `inbuilt_methods.jl`**: Manage configuration and assembly. The Builder translates flat, user-friendly configuration dictionaries into strictly typed `SimulationConfig` objects, while `inbuilt_methods.jl` defines pre-configured numerical schemes (e.g., Upwind, MUSCL with varying EPD MOOD strategies and a priori limiters).
*   **`SimulationFunctions/` & `time_integration.jl`**: Form the core execution engine. These modules orchestrate the main simulation loops, geometric boundary condition application, and time marching via explicit or IMEX Runge-Kutta stages.
*   **`macros.jl`**: Exposes the user-facing execution API (`@sim`, `@plot`, and `@backup`), enabling rapid simulation dispatch and seamless fallback to dynamic plugin loading.
*   **`Examples/`**: Contains standard, ready-to-run benchmark configurations (such as the Sod Shock Tube or Double Mach Reflection) demonstrating the framework's capabilities.

---

## Core Workflow

Running a simulation in HyperCloud requires zero boilerplate. Users define a setup function returning a `SimulationConfig` and utilize the macro API to execute it:

1.  **Define the Physics:** Create a configuration function specifying the geometry, boundary conditions, initial states, and selected numerical schemes.
2.  **Execute & Cache:** Run `@sim my_config` in the REPL. HyperCloud will automatically construct the grids, execute the time integration, and securely cache the serialized data via PDEStudioCore.
3.  **Visualize:** Run `@plot my_config` to launch the interactive PDEStudio GUI, instantly mapping your multidimensional meshfree data onto visual primitives.
4.  **Archive:** Use `@backup my_config` to save a timestamped copy of your exact configuration file, ensuring absolute reproducibility for publications.