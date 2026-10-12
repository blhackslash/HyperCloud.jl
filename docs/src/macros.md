# Execution Macros

The macros are the main user interface. They resolve a configuration by name, run it, and
handle plugin loading and world-age issues for you.

## Global settings

| Variable | Meaning |
|:---------|:--------|
| `CONFIG_SEARCH_PATHS` | Folders searched recursively for `name.jl` (default `["configs", "experiments"]`). |
| `CONFIG_ROOT_OVERWRITE` | Project root for that search. Defaults to `PDEStudioCore.get_save_path()`. |
| `BONITO_LISTEN_URL`, `BONITO_LISTEN_PORT`, `BONITO_PROXY_URL` | Server settings for browser plots with `WGLMakie`. |

All of them are `Ref`s, e.g. `HyperCloud.CONFIG_SEARCH_PATHS[] = ["my_cases"]`.

## `@sim`

```julia
@sim name
@sim "name"
@sim :name
```

Runs the configuration returned by `name()`:

1. If a function `name` exists in the calling module, it is used.
2. Otherwise the macro searches `CONFIG_SEARCH_PATHS` for `name.jl`, includes it, and calls the
   function of the same name.
3. The call goes through `Base.invokelatest`, so functions defined by a freshly included file
   work in the same session.

If neither exists, `@sim` throws an `ErrorException`.

!!! note "Failed runs do not abort a batch"
    The simulation runners catch exceptions, log them with `@error`, and return `nothing`.
    A batch therefore continues past a failing method; check the log for errors.

## `@plot`

```julia
@plot
@plot name
```

Opens the PDEStudio viewer. With `name`, the configuration is resolved like in `@sim`, passed
to the viewer with `set_sim_config!`, and computed (or loaded from cache).

The backend is chosen from what is loaded in the calling module:

- `GLMakie`: a local window.
- `WGLMakie`: plots are served by a Bonito server (configured once per session with the
  `BONITO_*` settings). If Bonito itself is not loaded, a warning is printed.
- Neither: a warning is printed and the figure is returned without being displayed.

## `@backup`

```julia
@backup name
```

Copies the file that defines `name` into a `backups` folder next to it, with a timestamp in the
file name (`name_YYYYMMDD_HHMMSS.jl`). Use it to archive the exact configuration behind a
result.

## API

```@docs
@sim
@plot
@backup
```
