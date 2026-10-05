export burgers_kinetic

function burgers_kinetic()
    varied_params = create_varied_dict()
    push!(varied_params, :Grid_SEED => [10, 100, 1000, 10000]) 
    
    base_methods = get_master_method_dict()
    
    # Generate symmetric kinetic velocities and indices for 1D, 1-component Burgers
    # v_mag = 2.0 (Must be strictly greater than max|u| of the IC to remain stable)
    kin_vels, kin_idx = generate_relax_params(2.0, 1, 1)
    
    sim_config = SimulationConfig(
        Dict{Symbol, Any}(
            # Global
            :snapshots => 25, 
            :remove_ghosts => false,
            :tmax => 2.0, # Shorter tmax for Burgers before the shock gets too steep
            
            # Time Namespace
            :Time_CFL => 0.2, 
            
            # Grid & Domain Namespace
            :Grid_domain => :rectangular,
            :Grid_Ns => (300,), 
            :Grid_mins => (-5.0,), 
            :Grid_maxs => (5.0,),
            :Grid_periodic => true,
            :Grid_randomness_factor => (0.2,),
            :Grid_seed => 42,
            :Grid_mover => :none,
            
            # Weight Namespace
            :Weight_name => :exponential,
            :Weight_alpha => 1.0,
            :Weight_range => 5.5,
            
            # IC Namespace (Gauss pulse of amplitude 1.0)
            :IC_name => :gauss,
            :IC_a => (1.0,),
            :IC_b => (0.0,),
            :IC_width => 1.0,
            
            # PDE Namespace
            :PDE_name => :burgers, 
            :PDE_D => 1,
            
            # Kinetic Namespace
            :Kinetic_velocities => kin_vels,
            :Kinetic_indices => kin_idx,
            :Kinetic_epsilon => 1e-6,

            :Scheme_MLS_order => 2,
        ),
        base_methods,
        # Compare a Direct scheme vs a Kinetic scheme
        ["RK2MUSCL2", "ARS222MUSCL2"];
        varied_params = varied_params, ref_func_name = "analytical_solution"
    )
    
    return sim_config
end