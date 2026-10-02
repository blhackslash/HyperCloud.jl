export advection_1d
function advection_1d()
    varied_params = create_varied_dict()
    # Ensure parameter variations use the updated flattened namespace keys
    push!(varied_params, :Grid_SEED => [10, 100, 1000, 10000]) 
    
    base_methods = get_master_method_dict()
    
    sim_config = SimulationConfig(
        Dict{Symbol, Any}(
            # Global
            :snapshots => 25, 
            :remove_ghosts => false,
            :sim_func_name => :run_simulation,
            :tmax => 10.0, 

            # Time Namespace
            :Time_CFL => 0.2, 
            
            # Grid & Domain Namespace
            :Grid_domain => :rectangular,
            :Grid_Ns => (1000,), 
            :Grid_mins => (-5.0,), 
            :Grid_maxs => (5.0,),
            :Grid_periodic => true,
            :Grid_randomness_factor => (0.2,),
            :Grid_SEED => 42,
            :Grid_mover => :none,
            
            # Weight Namespace
            :Weight_name => :exponential,
            :Weight_alpha => 1.0,
            :Weight_range => 5.5,
            
            # IC Namespace
            :IC_name => :gauss,
            :IC_a => (1.0,),
            :IC_b => (0.0,),
            :IC_width => 1.0,
            
            # PDE Namespace
            :PDE_name => :linear, 
            :PDE_velocities => ((1.0,),), 

            :Scheme_MLS_order => 2,
        ),
        base_methods,
        ["RK2MUSCL2"];
        varied_params = varied_params, ref_func_name = "analytical_solution"
    )
    return sim_config
end