export create_method, generate_relax_params, get_master_method_dict

# 1. Define the parameters that direct solvers should completely ignore
# Note: Updated to map to your new flattened namespace keys!
const ignore_relax = [:Kinetic_velocities, :Kinetic_epsilon, :Kinetic_indices, :Kinetic_interior_factor, :save_relax]

"""
    create_method(timestepper::Symbol, main_gradient::Symbol, order::Int; kwargs...)

Automatically generates a method ParamDict using the strict namespace API, 
automatically assigns the correct simulation pipeline, and generates a standardized plot label.
"""
function create_method(timestepper::Symbol, main_gradient::Symbol, order::Int;
                       main_flux::Symbol = :Rusanov,
                       limiter::Symbol = :none,
                       mood_criterion::Symbol = :none,
                       mood_strategy::Symbol = :EPD1,
                       sim_func_name::Union{Symbol, Nothing} = nothing,
                       custom_label::Union{String, Nothing} = nothing,
                       kwargs...)
                       
    # Automatically infer the simulation pipeline as a Symbol
    if isnothing(sim_func_name)
        is_imex = timestepper ∈ (:ARS222, :ARS233, :IMEXEuler, :PRSSP3, :SSP332)
        sim_func_name = is_imex ? :run_kinetic_simulation : :run_direct_simulation
    end

    # Build the base ParamDict mapped directly to the new API Namespaces
    dict = Dict{Symbol, Any}(
        :Time_stepper => timestepper,
        :Scheme_name => main_gradient,
        :Scheme_order => order,
        :Flux_name => main_flux,
        :Limiter_name => limiter,
        :MOOD_criterion => mood_criterion,
        :MOOD_strategy => mood_strategy,
        :sim_func_name => sim_func_name
    )
    
    # Add any extra kwargs (like the ignore lists)
    for (k, v) in kwargs
        dict[Symbol(k)] = v isa String ? Symbol(v) : v
    end
    
    # Auto-Generate the Label
    if !isnothing(custom_label)
        name = custom_label
    else
        name = "$(timestepper)$(main_gradient)$(order)"
        
        mods = String[]
        if main_flux != :Rusanov; push!(mods, string(main_flux)); end
        if limiter != :none; push!(mods, string(limiter)); end
        
        if mood_criterion != :none
            push!(mods, string(mood_strategy))
            if mood_criterion != :U2
                push!(mods, string(mood_criterion))
            end
        end
        
        if !isempty(mods)
            name *= "(" * join(mods, ",") * ")"
        end
    end
    
    return name => dict
end

"""
    get_master_method_dict()

Returns a MethodDict containing every standard numerical scheme in the Meshfree framework.
"""
function get_master_method_dict()
    methods = [
        "Analytical Solution" => Dict{Symbol, Any}(:sim_func_name => :run_direct_simulation),
        
        # =====================================================================
        # 1. STABLE BASELINES (1st Order)
        # =====================================================================
        create_method(:Euler, :Upwind, 1, ignore=ignore_relax),
        create_method(:RK2, :Upwind, 1, ignore=ignore_relax),
        create_method(:IMEXEuler, :Upwind, 1), 
        create_method(:ARS222, :Upwind, 1),
        create_method(:ARS233, :Upwind, 1),
        create_method(:RK2, :WENO, 2, ignore=ignore_relax),

        # =====================================================================
        # 2. EXPLICIT 2ND ORDER (RK2) - DIRECT
        # =====================================================================
        create_method(:RK2, :Upwind, 2, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  1, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, limiter=:VK, limiter_mode=:hard, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, limiter=:BJ, limiter_mode=:hard, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, limiter=:minmod, limiter_mode=:hard, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, mood_criterion=:U2, mood_strategy=:EPD0, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, mood_criterion=:U2, mood_strategy=:EPD1, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, mood_criterion=:U2, mood_strategy=:EPD2, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, mood_criterion=:U2, mood_strategy=:SEPD0, ignore=ignore_relax),
        create_method(:RK2, :MUSCL,  2, mood_criterion=:U1, ignore=ignore_relax),
        create_method(:RK4, :MUSCL,  4, mood_criterion=:only, ignore=ignore_relax),

        # =====================================================================
        # 3. RELAXATION 2ND ORDER (ARS222) - KINETIC
        # =====================================================================
        create_method(:ARS222, :Upwind, 2),
        create_method(:ARS222, :MUSCL,  2),
        create_method(:ARS222, :MUSCL,  2, limiter=:VK, limiter_mode=:hard),
        create_method(:ARS222, :MUSCL,  2, limiter=:minmod, limiter_mode=:hard),
        create_method(:ARS222, :MUSCL,  2, mood_criterion=:U2, mood_strategy=:EPD1),
        create_method(:ARS222, :MUSCL,  2, mood_criterion=:U2, mood_strategy=:EPD2),
        create_method(:ARS222, :MUSCL,  2, mood_criterion=:U2, mood_strategy=:SEPD0),
        create_method(:ARS222, :MUSCL,  2, limiter=:VK, limiter_mode=:hard),
        create_method(:ARS222, :MUSCL,  2, mood_criterion=:U1, mood_strategy=:EPD1),

        # =====================================================================
        # 4. HIGH-ORDER (RK4 Direct / ARS233 Kinetic)
        # =====================================================================
        (create_method(:RK4, :MUSCL, order, ignore=ignore_relax) for order in 2:5)...,
        (create_method(:RK4, :Upwind, order, ignore=ignore_relax) for order in 2:5)...,
        (create_method(:RK4, :MUSCL, order, limiter=:VK, limiter_mode=:hard, ignore=ignore_relax) for order in 2:5)...,
        (create_method(:RK4, :MUSCL, order, mood_criterion=:U2, mood_strategy=:EPD1, ignore=ignore_relax) for order in 2:5)...,
        (create_method(:RK4, :MUSCL, order, mood_criterion=:U2, mood_strategy=:EPD2, ignore=ignore_relax) for order in 2:5)...,
        (create_method(:ARS233, :MUSCL, order, mood_criterion=:U2, mood_strategy=:EPD1) for order in 2:5)...,
        (create_method(:RK4, :MUSCL, order, limiter=:minmod, limiter_mode=:hard, ignore=ignore_relax) for order in 2:5)...
    ]
    
    return Dict(methods...)
end

"""
    generate_relax_params(v_mag::Float64, dim::Int, NM::Int)

Generates symmetric kinetic velocities matching the native `Flux{D, NK}` tuple format 
along with the uniformly partitioned edge indices (`relax_indices`) for the kinetic system.
"""
function generate_relax_params(v_mag::Float64, dim::Int, NM::Int)
    kin_per_macro = 2 * dim
    NK = kin_per_macro * NM
    
    relax_velocities = ntuple(Val(dim)) do d
        ntuple(Val(NK)) do k
            m_idx = div(k - 1, kin_per_macro)
            local_k = k - m_idx * kin_per_macro
            
            d_active = div(local_k - 1, 2) + 1
            is_positive = isodd(local_k)
            
            if d == d_active
                return is_positive ? v_mag : -v_mag
            else
                return 0.0
            end
        end
    end
    
    relax_indices = zeros(Int, NM + 1)
    relax_indices[1] = 1
    for m in 1:NM
        relax_indices[m+1] = relax_indices[m] + kin_per_macro
    end
    
    return relax_velocities, relax_indices
end