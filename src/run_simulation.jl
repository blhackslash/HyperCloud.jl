export run_simulation


# ==============================================================================
# EXECUTION BARRIERS
# ==============================================================================
"""
    extract_namespace(params::Dict, prefix::Symbol; separator::String="_")

Extracts keys starting with `prefix` followed by `separator`. 
Leaves the values exactly as they are (preserving tuples for hashing).
Example: `:PDE_velocities` becomes `:velocities`.
"""
function extract_namespace(params::Dict, prefix::Symbol; separator::String="_")
    prefix_str = string(prefix) * separator
    len = length(prefix_str)
    
    extracted = Dict{Symbol, Any}()
    
    for (k, v) in params
        k_str = string(k)
        if startswith(k_str, prefix_str)
            # Slice off the prefix and keep the rest of the parameter name
            new_key = Symbol(k_str[len+1:end])
            extracted[new_key] = v
        end
    end
    
    return extracted
end

@inline _unwrap(v::SVector{1, T}) where {T} = v[1]
@inline _unwrap(v) = v

@noinline function _execute_explicit_sim!(method, eq, pg, dt, is_cfl, run_params, dimension, snapshots, remove_ghosts, M_components, xmins, xmaxs, ::Type{T}) where {T}
    xs_svector, us_svector, ts_full, k_step, elapsed_time = solve_equation(method, eq, pg, run_params[:tmax], dt; is_cfl = is_cfl, snapshots = snapshots, remove_ghosts = remove_ghosts, show_progress = _PROGRESS_BAR[], progress_interval = _PROGRESS_INTERVAL[])
    @info "Explicit Simulation (D=$dimension) finished in $(round(elapsed_time, digits=2)) seconds."

    valid_indices = findall(i -> isassigned(us_svector, i), 1:length(us_svector))
    ts = ts_full[valid_indices]
    
    xs_final = xs_svector[valid_indices]
    us_final = us_svector[valid_indices]
    tmin, tmax = zero(T), T(run_params[:tmax])
    
    sim_data_result = create_sim_data(xs_final, us_final, ts, run_params; xmins = xmins, xmaxs = xmaxs, tmin = tmin, tmax = tmax)
    add_stat!(sim_data_result, :Runtime, elapsed_time, :none)
    add_stat!(sim_data_result, :k_step, Float64(k_step), :none) 
    return sim_data_result
end

@noinline function _execute_kinetic_sim!(system_method, eq_kin, pg, dt, is_cfl, run_params, dimension, snapshots, remove_ghosts, save_relax, km, xmins, xmaxs, ::Type{T}) where {T}
    xs_svector, us_svector, ts_full, k_step, elapsed_time = solve_equation(system_method, eq_kin, pg, run_params[:tmax], dt; is_cfl = is_cfl, snapshots = snapshots, remove_ghosts = remove_ghosts, show_progress = _PROGRESS_BAR[], progress_interval = _PROGRESS_INTERVAL[])
    @info "Kinetic Relaxation Simulation (D=$dimension) finished in $(round(elapsed_time, digits=2)) seconds."

    valid_indices = findall(i -> isassigned(us_svector, i), 1:length(us_svector))
    ts = ts_full[valid_indices]
    xs_final = xs_svector[valid_indices]

    if save_relax
        us_final = us_svector[valid_indices]
    else
        us_final = [[km(u) for u in snap] for snap in us_svector[valid_indices]]
    end
    tmin, tmax = zero(T), T(run_params[:tmax])
    
    sim_data_result = create_sim_data(xs_final, us_final, ts, run_params; xmins = xmins, xmaxs = xmaxs, tmin = tmin, tmax = tmax)
    add_stat!(sim_data_result, :Runtime, elapsed_time, :none)
    add_stat!(sim_data_result, :k_step, Float64(k_step), :none) 
    return sim_data_result
end

# =========================================================================
# MODULAR TIME STEPPER BUILDERS
# =========================================================================

# --- Tableau Builder ---
function build_tableau(time_conf::Dict, context::Dict)
    name = time_conf[:stepper]::Symbol
    return build_tableau(Val(name), time_conf, context)
end

build_tableau(name::Val, conf::Dict, ctx::Dict) = error("Unknown Time Stepper Tableau: $(typeof(name))")

# Explicit RK Tableaus
build_tableau(::Val{:Euler}, conf::Dict, ctx::Dict) = RK1_Euler_Tableau(ctx[:Type])
build_tableau(::Val{:RK2}, conf::Dict, ctx::Dict)   = RK2_Ralston_Tableau(ctx[:Type])
build_tableau(::Val{:RK3}, conf::Dict, ctx::Dict)   = RK3_SSP_Tableau(ctx[:Type])
build_tableau(::Val{:RK4}, conf::Dict, ctx::Dict)   = RK4_Classical_Tableau(ctx[:Type])

# IMEX Tableaus
build_tableau(::Val{:ARS233}, conf::Dict, ctx::Dict) = IMEX_ARS233_Tableau(ctx[:Type])
build_tableau(::Val{:PRSSP3}, conf::Dict, ctx::Dict) = IMEX_PRSSP3_Tableau(ctx[:Type])
build_tableau(::Val{:ARS222}, conf::Dict, ctx::Dict) = IMEX_ARS222_Tableau(ctx[:Type])
build_tableau(::Val{:SSP332}, conf::Dict, ctx::Dict) = IMEX_SSP2332_Tableau(ctx[:Type])
build_tableau(::Val{:IMEXEuler}, conf::Dict, ctx::Dict) = IMEX_Euler_Tableau(ctx[:Type])

# --- Main Time Stepper Builder ---
function build_timestepper(context::Dict)
    tableau = context[:Tableau]
    return _build_timestepper(tableau, context)
end

# Deduce Explicit RK Stepper from the RKButcherTableau
function _build_timestepper(tableau::RKButcherTableau, context::Dict)
    return GeneralRKTimeStepper(
        context[:Equation], 
        context[:Scheme], 
        context[:ExplicitSources], 
        tableau
    )
end

# Deduce IMEX Stepper from the IMEXButcherTableau
function _build_timestepper(tableau::IMEXButcherTableau, context::Dict)
    return GeneralIMEXTimeStepper(
        context[:Equation], 
        context[:Scheme], 
        context[:ExplicitSources], 
        context[:ImplicitSources], 
        tableau
    )
end

# =========================================================================
# KINETIC SYSTEM BUILDER
# =========================================================================

function build_kinetic_system(kin_conf::Dict, context::Dict)
    T = context[:Type]::DataType
    D = context[:D]::Int
    NM = context[:M]::Int
    eq_macro = context[:Equation]
    
    # Check if kinetic relaxation is active
    if !haskey(kin_conf, :velocities)
        km = Kin2Macro(collect(1:(NM + 1)))
        return km, NM, nothing, nothing
    end
    
    relax_config = kin_conf[:velocities]::Tuple
    relax_eps = T(kin_conf[:epsilon])
    relax_indices = kin_conf[:indices]::Vector{Int}
    
    km = Kin2Macro(relax_indices)
    
    # eq_kin uses LinearAdvection (which safely triggers param2vel internally!)
    eq_kin = LinearAdvection(relax_config, T)
    
    # Extract exact number of kinetic components from the constructed PDE
    NK = typeof(eq_kin).parameters[2] 
    
    coeffs = State{NM, T}(ntuple(m -> one(T) / T(relax_indices[m+1] - relax_indices[m]), Val(NM)))
    interior_factor = T(get(kin_conf, :interior_factor, D)) 
    
    source_term = eq_macro.rep isa Conservative ? 
        RelaxationSourceTerm(km, relax_eps, coeffs, eq_macro, eq_kin, interior_factor) :
        NonLocalRelaxationSourceTerm(km, relax_eps, coeffs, eq_macro, eq_kin, interior_factor)
        
    return km, NK, eq_kin, source_term
end

function build_spatial_schemes(params::ParamDict, D::Int, M_comps::Int, delta_relax::T, ::Type{T}) where {T}
    order     = params[:order]
    main_grad = string(params[:main_gradient])
    main_flux = string(params[:main_flux])
    
    lim_name  = string(get(params, :limiter, "none"))
    lim_mode  = get(params, :limiter_mode, :soft)
    mood_crit = string(get(params, :mood_criterion, "none"))
    mood_strat= string(get(params, :mood_strategy, "EPD1"))
    mls_order = get(params, :MLS_order, 0)
    
    limiter = if lim_name == "minmod"; MinmodLimiter(lim_mode)
              elseif lim_name == "superbee"; SuperbeeLimiter(lim_mode)
              elseif lim_name == "VK"; VenkatakrishnanLimiter(lim_mode)
              elseif lim_name == "BJ"; BarthJespersenLimiter(lim_mode)
              else; NoLimiter() end

    mood_criterion = if mood_crit == "U2"; MOODu2(delta_relax)
               elseif mood_crit == "U1"; MOODu1(delta_relax)
               elseif mood_crit == "only"; OnlyMOOD()
               else; NoMOOD() end
               
    mood_strategy = if mood_strat == "EPD0"; EPD0()
                    elseif mood_strat == "SEPD0"; StrictEPD0()
                    elseif mood_strat == "EPD1"; EPD1()
                    elseif mood_strat == "EPD2"; EPD2()
                    end
    mood_fun = MOOD(mood_strategy, mood_criterion) 

    MainFlux = main_flux == "Rusanov" ? RusanovFlux() : (main_flux == "Upwind" ? UpwindFlux() : error("Flux NYI"))
    MainGrad = if main_grad == "MUSCL"
        MUSCL(T, D, M_comps, order; div_order = mls_order, flux = MainFlux, limiter = limiter, mood = mood_fun)
    elseif main_grad == "Upwind"
        UpwindDivergence(T, D, M_comps, order; flux=MainFlux, algType=get(params, :upwind_alg_nd, "Classic"))
    elseif main_grad == "Central"
        CentralDivergence(T, D, M_comps, order)
    elseif main_grad == "WENO"
        WENO(T, D, M_comps, order)
    else
        error("Main Gradient NYI") 
    end

    return MainGrad
end

# =========================================================================
# MODULAR SPATIAL SCHEME BUILDERS
# =========================================================================

# ==============================================================================
# MAIN SIMULATION ORCHESTRATOR
# ==============================================================================
function run_simulation(params::ParamDict)::Union{AbstractSimData, Nothing}
    @info "--- Running General N-Dimensional Simulation ---"
    
    # Base extraction mapping directly to Float64 by default
    T = get(params, :real_type, Float64) 
    
    try
        # Initialize the global build context
        context = Dict{Symbol, Any}()
        context[:Type] = T
        
        # 1. Extract Config Namespaces (Strictly)
        pde_conf    = extract_namespace(params, :PDE)
        domain_conf = extract_namespace(params, :Grid) 
        grid_conf   = extract_namespace(params, :Grid)
        weight_conf = extract_namespace(params, :Weight)
        flux_conf   = extract_namespace(params, :Flux)
        lim_conf    = extract_namespace(params, :Limiter)
        mood_conf   = extract_namespace(params, :MOOD)
        scheme_conf = extract_namespace(params, :Scheme)
        ic_conf     = extract_namespace(params, :IC)
        time_conf   = extract_namespace(params, :Time)
        
        # 2. Base Equation & Dimensions
        eq_macro = build_equation(pde_conf, context)
        
        D = typeof(eq_macro).parameters[1]
        NM = typeof(eq_macro).parameters[2]
        
        context[:D] = D
        context[:M] = NM
        
        # vel_var logic for grid movers
        vel_var = eq_macro isa EulerEquation ? Tuple(2:D+1) : (1,)
        
        # 3. Kinetic System Extension
        #km, M_comps, eq_kin, source_term = build_kinetic_system(params, D, NM, eq_macro.rep, T)
        
        # Update macro components to kinetic components if relaxation is active
        #context[:M] = M_comps
        #is_kinetic = !isnothing(eq_kin)
        is_kinetic = false
        
        # 4. Geometry & Particle Grid Configuration
        geom = build_domain(domain_conf, context)
        context[:Domain] = geom
        context[:WeightConf] = weight_conf
        
        pg = build_particle_grid(grid_conf, context)
        context[:Grid] = pg
        
        # 5. Pipeline State Registration
        context[:Equation] = is_kinetic ? eq_kin : eq_macro
        context[:ExplicitSources] = () # Empty tuple base for explicit sources
        #context[:ImplicitSources] = is_kinetic ? (source_term,) : ()
        
        # 6. Spatial Scheme Pipeline
        context[:Flux]    = build_flux(flux_conf, context)
        context[:Limiter] = build_limiter(lim_conf, context)
        context[:MOOD]    = build_mood(mood_conf, context)
        context[:Scheme]  = build_scheme(scheme_conf, context)
        
        # 7. Time Stepper Pipeline
        context[:Tableau] = build_tableau(time_conf, context)
        method = build_timestepper(context)
        
        # 8. Initial Condition Setup
        ic_conf = extract_namespace(params, :IC)
        IC = build_ic(ic_conf, context)
        if is_kinetic
            setInitialConditions!(pg, source_term, IC, eq_macro)
        else
            setInitialConditions!(pg, eq_macro, IC)
        end
        
        # 9. Time Constraints & Limits
        is_cfl = haskey(time_conf, :CFL)
        dt = is_cfl ? T(time_conf[:CFL]) : T(time_conf[:dt])
        
        geom_mins = Tuple(geom.mins)
        geom_maxs = Tuple(geom.maxs)
        
        # 10. Time Integration & Execution
        if !is_kinetic
            return _execute_explicit_sim!(method, eq_macro, pg, dt, is_cfl, params, D, params[:snapshots], get(params, :remove_ghosts, true), NM, geom_mins, geom_maxs, T)
        else
            return _execute_kinetic_sim!(method, eq_kin, pg, dt, is_cfl, params, D, params[:snapshots], get(params, :remove_ghosts, true), get(params, :save_relax, false), km, geom_mins, geom_maxs, T)
        end

    catch e
        @error "Error during Simulation!" exception=(e, catch_backtrace())
        return nothing
    end
end

function build_geometric_domain(params::ParamDict, D::Int, ::Type{T}) where {T}
    # 1. Parse Parameters & Catch Serialized Strings
    raw_bc = get(params, :bc, Dict{Int, AbstractBoundaryCondition}())
    
    if raw_bc isa String
        raw_bc = eval(Meta.parse(raw_bc))
    end
    
    # 2. Universal Boundary Condition Translation
    bc_map = Dict{Int, AbstractBoundaryCondition}()
    
    for (tag, bc_obj) in raw_bc
        if bc_obj isa AbstractBoundaryCondition
            bc_map[tag] = bc_obj
        else
            bc_sym = Symbol(bc_obj)
            if bc_sym === :outflow || bc_sym === :OutflowBC
                bc_map[tag] = OutflowBC()
            elseif bc_sym === :fixed_dirichlet || bc_sym === :FixedDirichlet
                bc_map[tag] = FixedDirichlet()
            else
                if isdefined(Main, bc_sym)
                    bc_map[tag] = getfield(Main, bc_sym)()
                elseif isdefined(@__MODULE__, bc_sym)
                    bc_map[tag] = getfield(@__MODULE__, bc_sym)()
                else
                    error("Boundary Condition '$bc_sym' could not be found.")
                end
            end
        end
    end
    
    # 3. Extract Periodic Flags for Geometry
    is_per_input = get(params, :periodic, false)
    is_per_svec = isa(is_per_input, Bool) ? SVector{D, Bool}(ntuple(_ -> is_per_input, Val(D))) : SVector{D, Bool}(is_per_input)
    
    # 4. Geometry Resolution
    domain_input = get(params, :domain, "rectangular")
    
    if typeof(domain_input) <: String || typeof(domain_input) <: Symbol
        shape = lowercase(string(domain_input))
        req_mins = T.(params[:mins]::Tuple)
        req_maxs = T.(params[:maxs]::Tuple)
        
        if shape == "rectangular"
            return get_rectangular_domain(T, req_mins, req_maxs; bc_map = bc_map, is_periodic = is_per_svec)
        elseif shape == "spherical"
            center = ntuple(d -> (req_mins[d] + req_maxs[d]) / 2.0, Val(D))
            radius = (req_maxs[1] - req_mins[1]) / 2.0
            return get_spherical_domain(T, center, radius; bc_map = bc_map, is_periodic = is_per_svec)
        else
            error("Unknown built-in domain shape: $shape")
        end
    else
        # Return the directly provided GeometricDomain instance
        return domain_input
    end
end

# =========================================================================
# INTERNAL PARTICLE GRID GENERATOR
# =========================================================================

function build_particle_grid(grid_conf::Dict, context::Dict)
    # 1. Pull required dependencies from the context
    T = context[:Type]::DataType
    D = context[:D]::Int
    M = context[:M]::Int
    geom = context[:Domain]
    
    geom_mins = Tuple(geom.mins)
    geom_maxs = Tuple(geom.maxs)
    
    # 2. Dynamically extract bounds from the resolved geometry
    Ns = grid_conf[:Ns]::Tuple
    rf_tuple = grid_conf[:randomness_factor]::Tuple
    
    dxs = ntuple(d -> (T(geom_maxs[d]) - T(geom_mins[d])) / Int(Ns[d]), Val(D))
    max_dx = maximum(dxs)
    nominal_dx = dxs
    randomness = ntuple(d -> T(rf_tuple[d]) * dxs[d], Val(D))
    
    # 3. Save max_dx into the context so the weight builder (and others) can access it
    context[:max_dx] = max_dx
    
    # 4. Instantiate the MLS Weight Function using the unified API
    weight_conf = context[:WeightConf]::Dict
    weight_func = build_weights(weight_conf, context)
    
    # Strict SEED extraction
    rng = MersenneTwister(grid_conf[:SEED]::Int)

    # 5. Construct the Particle Grid
    pg = ParticleGrid(
        geom, nominal_dx, weight_func, M;
        randomness = randomness,
        rng = rng,
    )
    
    return pg
end