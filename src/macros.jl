# Global configuration references for plotting
const BONITO_LISTEN_URL = Ref("0.0.0.0")
const BONITO_LISTEN_PORT = Ref(9384)
const BONITO_PROXY_URL = Ref("")
const BONITO_IS_CONFIGURED = Ref(false)

# Global configuration references for paths
const CONFIG_SEARCH_PATHS = Ref{Vector{String}}(["configs", "experiments"])
const CONFIG_ROOT_OVERWRITE = Ref{String}("")

"""
    _launch_and_display_plotter(mod::Module)

Internal helper function to check backends in the calling module, configure 
Bonito (if loaded), and trigger the main plot display. Returns `true` if successful.
"""
function _launch_and_display_plotter(mod::Module)
    if isdefined(mod, :launch_plotter) && isdefined(mod, :set_sim_config!)
        if isempty(BONITO_PROXY_URL[]); BONITO_PROXY_URL[] = "http://localhost:$(BONITO_LISTEN_PORT[])" end
        
        # --- Dynamic Backend Activation & Server Setup ---
        if isdefined(mod, :WGLMakie)
            if isdefined(mod, :Bonito)
                if !BONITO_IS_CONFIGURED[]
                    mod.Bonito.configure_server!(
                        listen_url = BONITO_LISTEN_URL[],
                        listen_port = BONITO_LISTEN_PORT[],
                        proxy_url = BONITO_PROXY_URL[]
                    )
                    BONITO_IS_CONFIGURED[] = true
                end
            else
                @warn "WGLMakie is loaded, but Bonito is not! Default Bonito settings are used. Please run `using Bonito` to configure the web server."
            end
            @info "Bonito server listens on $(BONITO_LISTEN_URL[]):$(BONITO_LISTEN_PORT[]). Visit $(BONITO_PROXY_URL[])/browser-display to see the plots." 
            mod.WGLMakie.activate!()
            
        elseif isdefined(mod, :GLMakie)
            mod.GLMakie.activate!()
        else
            @warn "No Makie backend found! Please run `using GLMakie` or `using WGLMakie`."
        end

        # --- Launch and Display ---
        fig = Base.invokelatest(() -> mod.launch_plotter())
        return fig
    else
        @warn "Please run `using PDEStudio` first to enable plotting."
        return nothing
    end
end

"""
    _get_user_project_root(mod::Module)

Attempts to find the root directory. Resolves in order:
1. `CONFIG_ROOT_OVERWRITE` (if explicitly set)
2. The calling module's project root (`pkgdir`)
3. The current working directory (`pwd()`)
"""
function _get_user_project_root(mod::Module)
    if !isempty(CONFIG_ROOT_OVERWRITE[])
        return CONFIG_ROOT_OVERWRITE[]
    end
    
    root = pkgdir(mod)
    return isnothing(root) ? pwd() : root
end

"""
    resolve_config_path(mod::Module, name::AbstractString)

Recursively searches for the configuration file `name.jl` across the user-defined `CONFIG_SEARCH_PATHS`.
"""
function resolve_config_path(mod::Module, name::AbstractString)
    root = _get_user_project_root(mod)
    target_file = "$(name).jl"
    
    for folder in CONFIG_SEARCH_PATHS[]
        search_base = joinpath(root, folder)
        if !isdir(search_base)
            continue
        end
        
        # Completely recursive search
        for (current_dir, dirs, files) in walkdir(search_base)
            if target_file in files
                found_path = joinpath(current_dir, target_file)
                @info "Loading plugin config from: $found_path"
                return found_path
            end
        end
    end
    
    error("Simulation setup for '$name' not found recursively in any of the configured search paths:$(CONFIG_SEARCH_PATHS[])")
end

"""
    backup_config(mod::Module, original_path::AbstractString, base_name::AbstractString)

Copies the executed file to a `backups` subfolder within the primary configuration directory.
"""
function backup_config(mod::Module, original_path::AbstractString, base_name::AbstractString)
    if !isfile(original_path)
        error("Source file for '$base_name' not found:$original_path")
    end

    root = _get_user_project_root(mod)
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    
    # Default to the first configured search path for saving backups
    backup_folder = isempty(CONFIG_SEARCH_PATHS[]) ? "configs" : first(CONFIG_SEARCH_PATHS[])
    backup_dir = joinpath(root, backup_folder, "backups")
    mkpath(backup_dir)
    
    backup_name = "$(base_name)_$(timestamp).jl"
    backup_path = joinpath(backup_dir, backup_name)
    
    cp(original_path, backup_path, force=true)
    @info "Archived config file for plotting/reproducibility: $backup_name"
    
    return backup_path
end

# Fallback method to extract the path directly from a Function object
function backup_config(mod::Module, func::Function, base_name::AbstractString)
    mths = methods(func)
    if isempty(mths)
        error("Could not find any methods for function '$base_name'.")
    end
    original_path = String(first(mths).file)
    return backup_config(mod, original_path, base_name)
end
# =========================================================================
# MACRO UTILITIES
# =========================================================================

# Helper to normalize macro inputs (bare words, strings, or symbols) into a pure Symbol
_normalize_macro_name(name) = name isa QuoteNode ? name.value : Symbol(name)

"""
    @sim name

Evaluates the configuration and runs the simulation. Checks the local module scope first; 
if not found, dynamically loads the configuration from the user's configured search paths.
Accepts strings (`"name"`), bare words (`name`), or symbols (`:name`).
"""
macro sim(name)
    mod = __module__
    raw_sym = _normalize_macro_name(name)
    name_sym = QuoteNode(raw_sym)
    name_str = string(raw_sym)
    
    return quote
        local config
        local is_local = isdefined($mod,$name_sym)
        
        if is_local
            # 1. Execute the function defined in the user's local scope
            config = $(esc(raw_sym))()
        else
            # 2. Plugin Fallback: Include the file directly into the user's module.
            path = HyperCloud.resolve_config_path($mod,$name_str)
            config = Base.include($mod, path)
            
            # If the file defines and returns a function instead of a raw config, evaluate it.
            if config isa Function
                config = config()
            end
        end
        
        Base.invokelatest(() -> run_all_simulations(config; force_overwrite=false, calculate_stats=true))
    end
end

"""
    @backup name

Manually triggers a backup of the source file where the configuration function or script `name` is defined.
"""
macro backup(name)
    mod = __module__
    raw_sym = _normalize_macro_name(name)
    name_sym = QuoteNode(raw_sym)
    name_str = string(raw_sym)

    return quote
        local is_local = isdefined($mod,$name_sym)
        
        if is_local
            HyperCloud.backup_config($mod, $(esc(raw_sym)),$name_str)
        else
            path = HyperCloud.resolve_config_path($mod,$name_str)
            HyperCloud.backup_config($mod, path,$name_str)
        end
    end
end

"""
    @plot()

Launches the PDEStudio plotter utilizing the active module's backend.
"""
macro plot()
    mod = __module__
    return quote
        HyperCloud._launch_and_display_plotter($mod)
    end
end

"""
    @plot name

Launches the PDEStudio plotter, loads the configuration from the function or plugin `name`, and forces a plot update.
"""
macro plot(name)
    mod = __module__
    raw_sym = _normalize_macro_name(name)
    name_sym = QuoteNode(raw_sym)
    name_str = string(raw_sym)

    return quote
        res = HyperCloud._launch_and_display_plotter($mod)
        if !isnothing(res)
            # Display GLMakie figure before changing it
            if !isdefined($mod, :WGLMakie); display(res) end
            
            local config
            local is_local = isdefined($mod,$name_sym)
            
            if is_local
                config = $(esc(raw_sym))()
            else
                path = HyperCloud.resolve_config_path($mod,$name_str)
                config = Base.include($mod, path)
                if config isa Function
                    config = config()
                end
            end

            $mod.set_sim_config!(config)
            $mod.force_simulation()
        else
            @warn "Plotter returned nothing!" 
        end
        res
    end
end