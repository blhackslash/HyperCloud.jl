export euler_1d_sod

# =========================================================================
# EXACT CONSTANTS FOR 1D SOD SHOCK TUBE
# =========================================================================

# Primitive states: (rho, u, p)
# Left: rho = 1.0, u = 0.0, p = 1.0
# E = p / (gamma - 1) + 0.5 * rho * u^2 = 1.0 / 0.4 = 2.5
const SOD_LEFT  = State{3, Float64}(1.0, 0.0, 2.5)

# Right: rho = 0.125, u = 0.0, p = 0.1
# E = 0.1 / 0.4 = 0.25
const SOD_RIGHT = State{3, Float64}(0.125, 0.0, 0.25)


# =========================================================================
# MAIN SIMULATION SETUP
# =========================================================================

function euler_1d_sod()
    # Interface normal for 1D
    p0 = (0.0,)
    n  = (1.0,) 
    
    # Simple boundaries for 1D Sod shock tube
    bc_map = Dict{Int, AbstractBoundaryCondition}(
        1 => FixedDirichlet(), # Left Boundary
        2 => OutflowBC()       # Right Boundary
    )

    # 1. Flat Namespace Shared Configuration
    shared = Dict{Symbol, Any}(
        :PDE_name => :euler,
        :PDE_D => 1,
        :PDE_gamma => 1.4,

        :sim_func_name => :run_direct_simulation,                
        :tmax => 0.2,
        :snapshots => 51,
        
        :Grid_domain => :rectangular,
        :Grid_mins => (-0.5,),
        :Grid_maxs => (0.5,),
        :Grid_Ns => (330,),
        :Grid_periodic => false,
        :Grid_bc_map => bc_map,
        :Grid_randomness_factor => (0.2,),
        :Grid_seed => 10,
        :Grid_mover => :none,
        
        # New Weight Function API mappings
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 5.5,
        
        :IC_name => :riemann,
        :IC_uL => SOD_LEFT,
        :IC_uR => SOD_RIGHT,
        :IC_p0 => p0,
        :IC_n => n,

        :Time_CFL => 0.2,
    )

    # 2. Modernized Method Configurations
    methods = Dict{Symbol, Dict{Symbol, Any}}(
        :analytical_solution => Dict(),
        # 1st Order Baseline
        :RK2Upwind1 => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :Upwind, :Scheme_order => 1, :Scheme_MLS_order => 1, 
            :Scheme_upwind_alg_nd => :Classic, :Flux_name => :Upwind,
            :MOOD_criterion => :none
        ),
        
        # A Priori Limiter
        :RK2MUSCL2_VK => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :MUSCL, :Scheme_order => 2, :Scheme_MLS_order => 0, 
            :Flux_name => :Rusanov, :Limiter_name => :VK, :Limiter_mode => :hard,
            :MOOD_criterion => :none
        ),
        
        # U2 MOOD Variations (with various EPD Halo Strategies)
        :RK2MUSCL2_U2_EPD1 => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :MUSCL, :Scheme_order => 2, :Scheme_MLS_order => 0, 
            :Flux_name => :Rusanov, :Limiter_name => :none,
            :MOOD_criterion => :U2, :MOOD_strategy => :EPD1, :MOOD_delta_relax => 1e-4
        ),
        :RK2MUSCL2_U2_EPD2 => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :MUSCL, :Scheme_order => 2, :Scheme_MLS_order => 0, 
            :Flux_name => :Rusanov, :Limiter_name => :none,
            :MOOD_criterion => :U2, :MOOD_strategy => :EPD2, :MOOD_delta_relax => 1e-4
        ),
        
        # U1 Standard DMP Variation
        :RK2MUSCL2_U1_EPD1 => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :MUSCL, :Scheme_order => 2, :Scheme_MLS_order => 0, 
            :Flux_name => :Rusanov, :Limiter_name => :none,
            :MOOD_criterion => :U1, :MOOD_strategy => :EPD1, :MOOD_delta_relax => 1e-4
        ),
        
        # High-Order RK4 MUSCL Variations
        :RK4MUSCL4_U2_EPD1 => Dict(
            :Time_stepper => :RK4, 
            :Scheme_name => :MUSCL, :Scheme_order => 4, :Scheme_MLS_order => 0, 
            :Flux_name => :Rusanov, :Limiter_name => :none,
            :MOOD_criterion => :U2, :MOOD_strategy => :EPD1, :MOOD_delta_relax => 1e-4
        )
    )

    # Active methods to run by default
    active_methods = [
        :analytical_solution,
        :RK2Upwind1,
        :RK2MUSCL2_VK,
        :RK2MUSCL2_U2_EPD1,
        :RK4MUSCL4_U2_EPD1
    ]

    return SimulationConfig(shared, methods, active_methods; ref_func_name=:analytical_solution)
end