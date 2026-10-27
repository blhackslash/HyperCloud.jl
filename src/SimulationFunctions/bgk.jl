
export run_bgk_simulation

"""
    run_bgk_simulation(params::ParamDict) -> Union{AbstractSimData, Nothing}

Runs the discrete-velocity BGK model (`:PDE_name => :bgk`, see [`BGKEquation`](@ref)). The particle
state is the distribution at the velocity nodes; the `IC` namespace gives the macroscopic state
`(ρ, ρu..., E)` (initial conditions are built with `M = D + 2`, like an Euler run with the same `γ`), and
every particle starts in its discrete Maxwellian. The collision term is implicit for an IMEX tableau
([`BGKRelaxation`](@ref)), explicit for an explicit one ([`ExplicitBGKRelaxation`](@ref)) and dropped in
the collisionless limit (`:PDE_tau => Inf`). Schemes, fluxes, limiters, MOOD and the grid are configured
exactly as for the other equations.

The returned snapshots hold the macroscopic moments `(ρ, ρu..., E)`, directly comparable with Euler
runs; with `:save_relax => true` they keep the full distributions instead.
`run_direct_simulation` and `run_kinetic_simulation` forward BGK configurations here, so the method
dictionaries work unchanged.
"""
function run_bgk_simulation(params::ParamDict)::Union{AbstractSimData, Nothing}
    T = get(params, :T, Float64)

    try
        @info "--- Running BGK Simulation ---"
        context = Dict{Symbol, Any}()
        context[:T] = T
        context[:Storage] = storage_for(get_backend())   # a GPU backend runs on its device
        apply_memory_limit!(context)                      # set_memory_limit!: reserve the buffers once

        # 1. Extract Config Namespaces
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
        adapt_conf  = extract_namespace(params, :Adapt)

        # 2. Equation; the initial condition lives in the macroscopic (Euler) variables
        eq = build_equation(pde_conf, context)
        eq isa BGKEquation || error("run_bgk_simulation needs a BGK equation (:PDE_name => :bgk), got $(typeof(eq))")
        D = get_D(eq)

        context[:Equation] = euler_limit(eq)
        context[:D] = D
        context[:M] = D + 2
        IC = build_ic(ic_conf, context)

        context[:Equation] = eq
        context[:M] = get_M(eq)

        # 3. Collision term: implicit for IMEX tableaus, explicit otherwise
        context[:Tableau] = build_tableau(time_conf, context)
        context[:ExplicitSources], context[:ImplicitSources] = bgk_sources(eq, context[:Tableau])
        collision = _is_collisionless(eq.relaxation) ? "none (collisionless)" :
                    isempty(context[:ImplicitSources]) ? "explicit" : "implicit"
        @info "BGK: D = $D, $(num_velocities(eq)) velocities, $(get_M(eq)) equations, γ = $(round(eq.gamma; digits = 4)), " *
              "collision: $collision"

        # 4. Geometry & Grid
        geom = build_domain(domain_conf, context)
        context[:Domain] = geom
        context[:WeightConf] = weight_conf

        pg = build_run_grid(grid_conf, context)           # also sets context[:Grid]
        adaptivity = build_adaptivity(adapt_conf, context)

        # 5. Spatial Schemes
        context[:Flux]    = build_flux(flux_conf, context)
        context[:Limiter] = build_limiter(lim_conf, context)
        context[:MOOD]    = build_mood(mood_conf, context)
        context[:Scheme]  = build_scheme(scheme_conf, context)

        # 6. Time Stepper
        method = build_run_timestepper(context; time_conf = time_conf)
        _log_memory(context)

        # 7. Initialization (discrete Maxwellian of the macroscopic initial state)
        set_initial_conditions!(pg, eq, IC)

        # 8. Time Settings
        is_cfl = haskey(time_conf, :CFL)
        dt = is_cfl ? T(time_conf[:CFL]) : T(time_conf[:dt])
        tmax = T(get(time_conf, :tmax, params[:tmax]))

        # 9. Execution
        xs_svector, us_svector, ts_full, k_step, elapsed_time = solve_equation(
            method, eq, pg, tmax, dt;
            is_cfl = is_cfl, snapshots = params[:snapshots],
            remove_ghosts = get(params, :remove_ghosts, true),
            show_progress = _PROGRESS_BAR[], progress_interval = _PROGRESS_INTERVAL[],
            adaptivity = adaptivity
        )

        @info "BGK Simulation (D=$D) finished in $(round(elapsed_time, digits=2)) seconds."

        # 10. Data Packaging: macroscopic moments unless the distributions are requested
        valid_indices = findall(i -> isassigned(us_svector, i), 1:length(us_svector))
        us_final = if get(params, :save_relax, false)
            us_svector[valid_indices]
        else
            [[moments(eq, u) for u in snap] for snap in us_svector[valid_indices]]
        end

        sim_data = create_sim_data(xs_svector[valid_indices], us_final, ts_full[valid_indices], params; xmins=Tuple(geom.mins), xmaxs=Tuple(geom.maxs), tmin=zero(T), tmax=tmax)
        add_stat!(sim_data, :Runtime, elapsed_time, :none)
        add_stat!(sim_data, :k_step, Float64(k_step), :none)
        add_stat!(sim_data, :final_particle_count, Float64(pg.meta.N), :none)
        record_memory!(sim_data, context)

        return sim_data

    catch e
        @error "Error during BGK Simulation!" exception=(e, catch_backtrace())
        return nothing
    end
end

# BGK configurations carry their own collision term, so the direct and kinetic runners hand them over.
_is_bgk(params) = get(params, :PDE_name, nothing) === :bgk
