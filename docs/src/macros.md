# Execution Macros

The execution macros form the primary user interface of the HyperCloud framework, eliminating boilerplate code and seamlessly handling module scopes, dynamic plugin loading, and world-age resolution.

## Global Configuration Variables

The macro system is governed by several global `Ref` variables that dictate file search paths and plotting behavior:

*   `CONFIG_SEARCH_PATHS`: A vector of folders (default: `["configs", "experiments"]`) where the macro recursively searches for configuration scripts.
*   `CONFIG_ROOT_OVERWRITE`: An explicit override for the project root directory. The default is `PDEStudioCore.get_save_path()`.
*   `BONITO_LISTEN_URL`, `BONITO_LISTEN_PORT`, `BONITO_PROXY_URL`: Network configurations for hosting interactive WGLMakie plots via a headless Bonito server.

## Macro API

### `@sim`

```julia
@sim name
```

The core execution macro. It triggers the construction and integration of the simulation configuration.

1.  **Scope Resolution:** It first checks if the function `name` exists in the local module scope.
2.  **Plugin Fallback:** If not found locally, it recursively searches the user's `CONFIG_SEARCH_PATHS` for a file named `name.jl`, loads it dynamically via `Base.include`, and executes the returned `SimulationConfig`.
3.  **World Age Resolution:** It uses `Base.invokelatest` to bypass Julia's strict world-age limitations, ensuring dynamically loaded plugins compile and execute flawlessly.

### `@plot`

```julia
@plot()
@plot name
```

The plotting macro bridges the computational backend with the interactive `PDEStudio` visualization frontend.

*   **Backend Activation:** The `_launch_and_display_plotter` helper dynamically detects if `GLMakie` or `WGLMakie` is loaded in the calling module. If `WGLMakie` is active, it automatically configures a Bonito server to serve the interactive plots over the network.
*   **Execution:** When called with `name`, it resolves the configuration exactly like `@sim`, injects it into the PDEStudio UI via `set_sim_config!`, and forces an immediate visual update.

### `@backup`

```julia
@backup name
```

Ensures exact reproducibility. It locates the source file where `name` is defined and explicitly copies it into a timestamped `backups` subdirectory.

## Documentation

```@docs
@sim
@plot
@backup
```