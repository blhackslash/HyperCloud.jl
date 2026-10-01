# 1. Define the parameters that direct solvers should completely ignore
const ignore_relax = [:relax_velocities, :relax_epsilon, :save_relax, :interior_factor, :relax_indices]

"""
    create_method(timestepper::String, main_gradient::String, order::Int; kwargs...)

Automatically generates a method ParamDict using the strict namespace API and a standardized plot label[cite: 31].
"""
function create_method(timestepper::String, main_gradient::String, order::Int;
                       main_flux::String = "Rusanov",
                       limiter::String = "none",
                       mood_criterion::String = "none",
                       mood_strategy::String = "EPD1",
                       custom_label::Union{String, Nothing} = nothing,
                       kwargs...)
                       
    # 1. Build the base ParamDict mapped directly to the new API Namespaces
    dict = Dict{Symbol, Any}(
        :Time_stepper => Symbol(timestepper),
        :Scheme_name => Symbol(main_gradient),
        :Scheme_order => order,
        :Flux_name => Symbol(main_flux),
        :Limiter_name => Symbol(limiter),
        :MOOD_criterion => Symbol(mood_criterion),
        :MOOD_strategy => Symbol(mood_strategy),
    )
    
    # 2. Add any extra kwargs (ignores, PDE overrides, relax_epsilon, etc.)
    for (k, v) in kwargs
        # Automatically cast string values to Symbols if they are added manually
        dict[Symbol(k)] = v isa String ? Symbol(v) : v
    end
    
    # 3. Auto-Generate the Label (if a custom one wasn't forced)
    if !isnothing(custom_label)
        name = custom_label
    else
        name = "$(timestepper)$(main_gradient)$(order)"
        
        mods = String[]
        if main_flux != "Rusanov"; push!(mods, main_flux); end
        if limiter != "none"; push!(mods, limiter); end
        
        if mood_criterion != "none"
            push!(mods, mood_strategy)
            if mood_criterion != "U2"
                push!(mods, mood_criterion)
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

Returns a MethodDict containing every standard numerical scheme in the Meshfree framework[cite: 29].
Direct methods automatically ignore relaxation-specific shared parameters[cite: 29].
"""
function get_master_method_dict()
    methods = [
        "Analytical Solution" => Dict{Symbol, Any}(),
        
        # 1. STABLE BASELINES (1st Order)
        create_method("Euler", "Upwind", 1, ignore=ignore_relax),
        create_method("RK2", "Upwind", 1, ignore=ignore_relax),
        create_method("IMEXEuler", "Upwind", 1), 
        create_method("ARS222", "Upwind", 1),
        create_method("ARS233", "Upwind", 1),
        create_method("RK2", "WENO", 2, ignore=ignore_relax),

        # 2. EXPLICIT 2ND ORDER (RK2)
        create_method("RK2", "Upwind", 2, ignore=ignore_relax),
        create_method("RK2", "MUSCL",  2, ignore=ignore_relax),
        create_method("RK2", "MUSCL",  1, ignore=ignore_relax),
        create_method("RK2", "MUSCL",  2, limiter="VK", ignore=ignore_relax, limiter_mode=:hard),
        create_method("RK2", "MUSCL",  2, limiter="BJ", ignore=ignore_relax, limiter_mode=:hard),
        create_method("RK2", "MUSCL",  2, limiter="minmod", ignore=ignore_relax, limiter_mode=:hard),
        create_method("RK2", "MUSCL",  2, mood_criterion="U2", mood_strategy="EPD0", ignore=ignore_relax),
        create_method("RK2", "MUSCL",  2, mood_criterion="U2", mood_strategy="EPD1", ignore=ignore_relax),
        create_method("RK2", "MUSCL",  2, mood_criterion="U2", mood_strategy="EPD2", ignore=ignore_relax),
        create_method("RK2", "MUSCL",  2, mood_criterion="U2", mood_strategy="SEPD0", ignore=ignore_relax),
        create_method("RK2", "MUSCL",  2, mood_criterion="U1", ignore=ignore_relax),
        create_method("RK4", "MUSCL",  4, mood_criterion="only", ignore=ignore_relax),

        # 3. RELAXATION 2ND ORDER (ARS222)
        create_method("ARS222", "Upwind", 2),
        create_method("ARS222", "MUSCL",  2),
        create_method("ARS222", "MUSCL",  2, limiter="VK", limiter_mode=:hard),
        create_method("ARS222", "MUSCL",  2, limiter="minmod", limiter_mode=:hard),
        create_method("ARS222", "MUSCL",  2, mood_criterion="U2", mood_strategy="EPD1", ignore=ignore_relax),
        create_method("ARS222", "MUSCL",  2, mood_criterion="U2", mood_strategy="EPD2", ignore=ignore_relax),
        create_method("ARS222", "MUSCL",  2, mood_criterion="U2", mood_strategy="SEPD0", ignore=ignore_relax),
        create_method("ARS222", "MUSCL",  2, limiter="VK", ignore=ignore_relax, limiter_mode=:hard),
        create_method("ARS222", "MUSCL",  2, mood_criterion="U1", mood_strategy="EPD1"),

        # 4. HIGH-ORDER EXPLICIT DIRECT (RK4)
        (create_method("RK4", "MUSCL", order, ignore=ignore_relax) for order in 2:5)...,
        (create_method("RK4", "Upwind", order, ignore=ignore_relax) for order in 2:5)...,
        (create_method("RK4", "MUSCL", order, limiter="VK", ignore=ignore_relax, limiter_mode=:hard) for order in 2:5)...,
        (create_method("RK4", "MUSCL", order, mood_criterion="U2", mood_strategy="EPD1", ignore=ignore_relax) for order in 2:5)...,
        (create_method("RK4", "MUSCL", order, mood_criterion="U2", mood_strategy="EPD2", ignore=ignore_relax) for order in 2:5)...,
        (create_method("ARS233", "MUSCL", order, mood_criterion="U2", mood_strategy="EPD1") for order in 2:5)...,
        (create_method("RK4", "MUSCL", order, limiter="minmod", ignore=ignore_relax, limiter_mode=:hard) for order in 2:5)...
    ]
    
    return Dict(methods...)
end

"""
    generate_relax_params(v_mag::Float64, dim::Int, NM::Int)

Generates symmetric kinetic velocities matching the native `Flux{D, NK}` tuple format 
along with the uniformly partitioned edge indices (`relax_indices`) for the kinetic system[cite: 31].
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