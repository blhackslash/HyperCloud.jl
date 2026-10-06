

using HyperCloud
using StaticArrays
using Test
using LinearAlgebra
using Dates

# =========================================================================
# 0. MOCK MODULES (Must be defined at the top-level scope)
# =========================================================================

module MockPDEStudio
    using HyperCloud
    using Base
    
    # Mock the required PDEStudio hooks
    launch_plotter() = "mock_figure"
    set_sim_config!(conf) = nothing
    force_simulation() = nothing
    
    # Dummy config for the macro test
    dummy_config() = "mock_config"
end

module MockGLMakie
    activate!() = nothing
end

module MockWGLMakie
    activate!() = nothing
end

module MockBonito
    configure_server!(; listen_url, listen_port, proxy_url) = nothing
end

# =========================================================================
# MAIN TEST SUITE
# =========================================================================
@testset "HyperCloud Tests" begin
@testset "Macro Utilities & Simulation Launchers" begin
    # ---------------------------------------------------------
    # 1. SETUP ISOLATED TEST ENVIRONMENT
    # ---------------------------------------------------------
    test_dir = mktempdir()
    mock_configs_dir = joinpath(test_dir, "mock_configs")
    mkdir(mock_configs_dir)
    
    # Overwrite the root so the macros search our temp directory
    HyperCloud.CONFIG_ROOT_OVERWRITE[] = test_dir
    HyperCloud.CONFIG_SEARCH_PATHS[] = ["mock_configs"]
    
    # ---------------------------------------------------------
    # 2. MOCK CONFIGURATIONS & FUNCTIONS
    # ---------------------------------------------------------
    
    # A. A dummy local function to test in-scope @sim
    @eval Main begin
        function my_local_sim()
            shared = Dict{Symbol, Any}(:Time_tmax => 0.0, :Time_snapshots => 0)
            return SimulationConfig(shared, Dict(), Symbol[])
        end
    end
    
    # B. A mock plugin script written to the filesystem
    plugin_name = "test_plugin"
    plugin_path = joinpath(mock_configs_dir, "$(plugin_name).jl")
    open(plugin_path, "w") do f
        write(f, """
        function test_plugin()
            shared = Dict{Symbol, Any}(:Time_tmax => 0.0, :Time_snapshots => 0)
            return SimulationConfig(shared, Dict(), Symbol[])
        end
        """)
    end
    
    # ---------------------------------------------------------
    # 3. TEST @sim MACRO (Local vs Plugin Fallback)
    # ---------------------------------------------------------
    @testset "@sim Macro Execution" begin
        # 1. Should successfully run the local function defined above
        @eval Main @sim my_local_sim
        @eval Main @sim "my_local_sim"
        @eval Main @sim :my_local_sim
        @test true
        
        # 2. Should dynamically find and run the plugin we wrote to disk
        @eval Main @sim test_plugin
        @test true
        
        # 3. Should fail gracefully if neither exists
        @test_throws ErrorException @eval Main @sim missing_sim_config
    end
    
    # ---------------------------------------------------------
    # 4. TEST @backup MACRO
    # ---------------------------------------------------------
    @testset "@backup Macro" begin
        # Backup the dynamically loaded plugin
        @eval Main @backup test_plugin
        
        # Check that the backup folder and a timestamped file were created
        backup_dir = joinpath(mock_configs_dir, "backups")
        @test isdir(backup_dir)
        
        files = readdir(backup_dir)
        @test length(files) == 1
        @test startswith(files[1], "test_plugin_")
        @test endswith(files[1], ".jl")
    end
    
    # ---------------------------------------------------------
    # 5. TEST @plot MACROS & BACKEND ROUTING
    # ---------------------------------------------------------
    @testset "@plot Macro (Backend Routing)" begin
        # Reset Bonito configuration state
        HyperCloud.BONITO_IS_CONFIGURED[] = false
        
        # Test 1: No Backend Loaded 
        # (It throws a warning, but still returns the figure)
        res = @test_logs (:warn, r"No Makie backend found!") HyperCloud._launch_and_display_plotter(MockPDEStudio)
        @test res == "mock_figure"
        
        # Test 2: GLMakie Simulation
        # Inject the mock GLMakie into MockPDEStudio
        @eval MockPDEStudio GLMakie = $MockGLMakie
        res = HyperCloud._launch_and_display_plotter(MockPDEStudio)
        @test res == "mock_figure"
        
        # Test 3: WGLMakie Simulation (Without Bonito)
        @eval MockPDEStudio WGLMakie = $MockWGLMakie
        res = @test_logs (:warn, r"WGLMakie is loaded, but Bonito is not!") (:info, r"Bonito server listens on") HyperCloud._launch_and_display_plotter(MockPDEStudio)
        @test res == "mock_figure"
        
        # Test 4: WGLMakie + Bonito Configuration
        @eval MockPDEStudio Bonito = $MockBonito
        @test HyperCloud.BONITO_IS_CONFIGURED[] == false
        res = @test_logs (:info, r"Bonito server listens on") HyperCloud._launch_and_display_plotter(MockPDEStudio)
        @test HyperCloud.BONITO_IS_CONFIGURED[] == true
        
        # Test 5: Full @plot Macro Execution within the Mock Module
        res = @eval MockPDEStudio @plot dummy_config
        @test res == "mock_figure"
    end
end

@testset "Initial Conditions - 2D Linear Advection" begin
    # =========================================================================
    # SHARED SIMULATION CONFIGURATION
    # =========================================================================
    
    # Advect diagonally to test multi-dimensional boundary handling
    velocities = ((1.,),(1.,))
    
    # Common simulation parameters
    shared = Dict{Symbol, Any}(
        :PDE_name => :linear,
        :PDE_D => 2,
        :PDE_velocities => velocities,
        
        :sim_func_name => :run_direct_simulation,                
        :tmax => 0.1, # Short duration just to verify integration works
        :snapshots => 2,
        
        :Grid_domain => :rectangular,
        :Grid_mins => (0.0, 0.0),
        :Grid_maxs => (1.0, 1.0),
        :Grid_Ns => (20, 20),
        :Grid_periodic => true,
        :Grid_randomness_factor => (0.0, 0.0), # Perfect grid for analytic testing
        :Grid_seed => 42,
        :Grid_mover => :none,
        
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 2.5,

        :Time_CFL => 0.2,
    )
    
    # We will use the simplest Upwind method to drive the integration
    methods = Dict{Symbol, Dict{Symbol, Any}}(
        :Upwind => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :Upwind, :Scheme_order => 1, :Scheme_MLS_order => 1, 
            :Scheme_upwind_alg_nd => :Classic, :Flux_name => :Upwind,
            :MOOD_criterion => :none
        )
    )
    
# =========================================================================
    # TESTS
    # =========================================================================

    @testset "IC: Box" begin
        ic_dict = Dict(
            :IC_name => :box,
            :IC_u_bg => (0.1,),
            :IC_u_box => (1.0,),
            :IC_mins => (0.25, 0.25),
            :IC_maxs => (0.75, 0.75)
        )
        
        conf = SimulationConfig(merge(shared, ic_dict), methods, [:Upwind])
        run_all_simulations(conf; force_overwrite=true, calculate_stats=false)
        @test true 
    end

    @testset "IC: Gauss" begin
        ic_dict = Dict(
            :IC_name => :gauss,
            :IC_a => (1.0,),
            :IC_b => (0.5, 0.5),
            :IC_width => 0.1
        )
        
        conf = SimulationConfig(merge(shared, ic_dict), methods, [:Upwind])
        run_all_simulations(conf; force_overwrite=true, calculate_stats=false)
        @test true
    end

    @testset "IC: Sine" begin
        ic_dict = Dict(
            :IC_name => :sine,
            :IC_a => (1.0,),
            :IC_period => (1.0, 1.0),
            :IC_c_offset => (0.5,)
        )
        
        conf = SimulationConfig(merge(shared, ic_dict), methods, [:Upwind])
        run_all_simulations(conf; force_overwrite=true, calculate_stats=false)
        @test true
    end

    @testset "IC: Standard Riemann" begin
        ic_dict = Dict(
            :IC_name => :riemann,
            :IC_uL => (1.0,),
            :IC_uR => (0.0,),
            :IC_p0 => (0.5, 0.5),
            :IC_n => (1.0, 1.0) # Diagonal normal
        )
        
        # We need to turn off periodic boundaries for a non-periodic Riemann problem
        bc_map = Dict{Int, AbstractBoundaryCondition}(
            1 => OutflowBC(), 2 => OutflowBC(), 3 => OutflowBC(), 4 => OutflowBC()
        )
        local_shared = copy(shared)
        local_shared[:Grid_periodic] = false
        local_shared[:Grid_bc_map] = bc_map
        
        conf = SimulationConfig(merge(local_shared, ic_dict), methods, [:Upwind])
        run_all_simulations(conf; force_overwrite=true, calculate_stats=false)
        @test true
    end

    @testset "IC: Smoothed Riemann" begin
        ic_dict = Dict(
            :IC_name => :s_riemann,
            :IC_uL => (1.0,),
            :IC_uR => (0.0,),
            :IC_p0 => (0.5, 0.5),
            :IC_n => (1.0, 0.0),
            :IC_width => 0.05
        )
        
        bc_map = Dict{Int, AbstractBoundaryCondition}(
            1 => OutflowBC(), 2 => OutflowBC(), 3 => OutflowBC(), 4 => OutflowBC()
        )
        local_shared = copy(shared)
        local_shared[:Grid_periodic] = false
        local_shared[:Grid_bc_map] = bc_map
        
        conf = SimulationConfig(merge(local_shared, ic_dict), methods, [:Upwind])
        run_all_simulations(conf; force_overwrite=true, calculate_stats=false)
        @test true
    end

    @testset "IC: Quadrant Riemann (QRiemann)" begin
        states = (
            (1.0,), # Bottom-Left
            (2.0,), # Bottom-Right
            (3.0,), # Top-Left
            (4.0,)  # Top-Right
        )
        
        ic_dict = Dict(
            :IC_name => :q_riemann,
            :IC_states => states,
            :IC_p0 => (0.5, 0.5)
        )
        
        bc_map = Dict{Int, AbstractBoundaryCondition}(
            1 => OutflowBC(), 2 => OutflowBC(), 3 => OutflowBC(), 4 => OutflowBC()
        )
        local_shared = copy(shared)
        local_shared[:Grid_periodic] = false
        local_shared[:Grid_bc_map] = bc_map
        
        conf = SimulationConfig(merge(local_shared, ic_dict), methods, [:Upwind])
        run_all_simulations(conf; force_overwrite=true, calculate_stats=false)
        @test true
    end

end
@testset "Analytical Closures & Solutions" begin

    # Helper to quickly generate a real Rectangular Domain using your builder
    function get_test_domain(mins, maxs, periodic)
        D = length(mins)
        conf = Dict{Symbol, Any}(
            :domain => :rectangular, 
            :mins => mins, 
            :maxs => maxs, 
            :periodic => periodic
        )
        ctx = Dict{Symbol, Any}(:T => Float64, :D => D)
        return build_domain(conf, ctx)
    end

    @testset "Linear Advection (Characteristic Tracing)" begin
        # 1D Advection with velocity v = 2.0
        eq = LinearAdvection(((2.0,),), Float64)
        ic = Gauss(State{1, Float64}(1.0), Space{1, Float64}(0.0), 1.0)
        
        geom = get_test_domain((-5.0,), (5.0,), true)
        exact_func = analytic_closure(eq, ic, geom)
        
        # At t = 0, the analytical solution should perfectly match the IC
        st_t0 = SVector{2, Float64}(0.0, 0.0) # (x, t)
        @test exact_func(st_t0)[1] ≈ 1.0
        
        # At t = 1.0, the peak should have moved to x = 2.0
        st_t1_peak = SVector{2, Float64}(2.0, 1.0)
        @test exact_func(st_t1_peak)[1] ≈ 1.0
        
        # Periodicity Check: At t = 5.0, peak moved 10 units -> wrapped back to x = 0.0
        st_t5_peak = SVector{2, Float64}(0.0, 5.0)
        @test exact_func(st_t5_peak)[1] ≈ 1.0
    end

    @testset "Burgers Equation" begin
        eq = BurgersEquation{1, Float64}()
        geom = get_test_domain((-5.0,), (5.0,), false)
        
        @testset "Riemann Shock" begin
            # uL = 2.0, uR = 0.0 -> Shock speed s = 1.0
            ic_shock = Riemann(State{1}(2.0), State{1}(0.0), Space{1}(0.0), Space{1}(1.0))
            exact_shock = analytic_closure(eq, ic_shock, geom)
            
            # At t = 1.0, the shock is at x = 1.0
            @test exact_shock(SVector{2, Float64}(0.5, 1.0))[1] ≈ 2.0 # Left of shock
            @test exact_shock(SVector{2, Float64}(1.5, 1.0))[1] ≈ 0.0 # Right of shock
        end
        
        @testset "Riemann Rarefaction" begin
            # uL = 0.0, uR = 2.0 -> Rarefaction fan
            ic_rare = Riemann(State{1}(0.0), State{1}(2.0), Space{1}(0.0), Space{1}(1.0))
            exact_rare = analytic_closure(eq, ic_rare, geom)
            
            # At t = 1.0, fan spans from x = 0.0 to x = 2.0
            @test exact_rare(SVector{2, Float64}(-1.0, 1.0))[1] ≈ 0.0 # Left of fan
            @test exact_rare(SVector{2, Float64}(1.0, 1.0))[1]  ≈ 1.0 # Inside fan (x/t = 1.0/1.0)
            @test exact_rare(SVector{2, Float64}(3.0, 1.0))[1]  ≈ 2.0 # Right of fan
        end
    end

    @testset "Euler Equation (Sod Shock Tube)" begin
        eq = EulerEquation(Val(1), Float64, 1.4, Conservative())
        geom = get_test_domain((-5.0,), (5.0,), false)
        
        # Sod Shock Tube States
        SOD_LEFT  = State{3, Float64}(1.0, 0.0, 2.5)
        SOD_RIGHT = State{3, Float64}(0.125, 0.0, 0.25)
        ic_sod = Riemann(SOD_LEFT, SOD_RIGHT, Space{1}(0.0), Space{1}(1.0))
        
        exact_sod = analytic_closure(eq, ic_sod, geom)
        
        # At t = 0, check the unperturbed states
        @test exact_sod(SVector{2, Float64}(-1.0, 0.0))[1] ≈ 1.0
        @test exact_sod(SVector{2, Float64}(1.0, 0.0))[1] ≈ 0.125
        
        # At t = 0.2, the right shock should have moved rightwards (density > 0.125)
        state_mid = exact_sod(SVector{2, Float64}(0.1, 0.2))
        @test 0.125 < state_mid[1] < 1.0 # Density in the star region
    end

    @testset "Numerical vs Analytical Consistency" begin
        # 1. Grab the config generator for Sod Shock Tube
        config = euler_1d_sod()
        
        # 2. Extract the baseline explicit RK2 scheme dict and merge it with shared params
        # This gives us a complete flat ParamDict ready for execution
        baseline_method_name = :RK2Upwind1
        params = merge(config.shared_params, config.methods_dict[baseline_method_name])
        
        # 3. Build the exact analytical closure using the analytical_solution wrapper
        exact_u = analytical_solution(params)
        
        # 4. Execute the numerical simulation headlessly
        sim_data = run_direct_simulation(params)
        
        # 5. Extract the final snapshot
        final_xs = sim_data.x[end]
        final_us = sim_data.u[end]
        final_t  = sim_data.t[end]
        
        error_L1 = 0.0
        for i in eachindex(final_xs)
            st = SVector{2, Float64}(final_xs[i][1], final_t)
            u_exact = exact_u(st)
            u_num = final_us[i]
            
            # Sum absolute density error
            error_L1 += abs(u_exact[1] - u_num[1])
        end
        error_L1 /= length(final_xs)
        
        # 6. Evaluate if the average density error is within reasonable bounds
        # An L1 error of < 0.1 confirms that the primary numerical framework logic is sound
        @test error_L1 < 0.1 
    end
    @testset "2D Burgers Equation" begin
        # Instantiate the 2D Burgers Equation (M = 1 strictly enforced)
        eq_2d = BurgersEquation{2, Float64}()
        geom_2d = get_test_domain((-5.0, -5.0), (5.0, 5.0), (false, false))
        
        @testset "Gauss (Characteristic Trace)" begin
            # 2D Gauss IC: amplitude 1.0, centered at (0,0), width 1.0
            ic_gauss = Gauss(State{1}(1.0), Space{2}(0.0, 0.0), 1.0)
            exact_gauss = analytic_closure(eq_2d, ic_gauss, geom_2d)
            
            # At t = 0, check the center peak and a spatial offset
            @test exact_gauss(SVector{3, Float64}(0.0, 0.0, 0.0))[1] ≈ 1.0
            @test exact_gauss(SVector{3, Float64}(1.0, 0.0, 0.0))[1] ≈ exp(-1.0)
            
            # At t = 0.5, the peak (u=1) advects isotropically at velocity (1, 1)
            # Therefore, the peak should be exactly at (0.5, 0.5)
            @test exact_gauss(SVector{3, Float64}(0.5, 0.5, 0.5))[1] ≈ 1.0
        end
        
        @testset "Riemann Shock" begin
            # Diagonal shock wave
            n_diag = normalize(Space{2}(1.0, 1.0))
            ic_shock = Riemann(State{1}(2.0), State{1}(0.0), Space{2}(0.0, 0.0), n_diag)
            exact_shock = analytic_closure(eq_2d, ic_shock, geom_2d)
            
            # n_sum = sqrt(2). uL = 2.0, uR = 0.0 -> planar shock speed s = sqrt(2)
            # At t = 1.0, the planar shock front is at distance d = sqrt(2) along the normal
            
            # Test a point strictly behind the shock front (d < sqrt(2))
            @test exact_shock(SVector{3, Float64}(0.5, 0.5, 1.0))[1] ≈ 2.0
            
            # Test a point strictly ahead of the shock front (d > sqrt(2))
            @test exact_shock(SVector{3, Float64}(2.0, 2.0, 1.0))[1] ≈ 0.0
        end

        @testset "Riemann Rarefaction" begin
            # 2D Rarefaction along the X axis
            n_x = Space{2}(1.0, 0.0)
            ic_rare = Riemann(State{1}(0.0), State{1}(2.0), Space{2}(0.0, 0.0), n_x)
            exact_rare = analytic_closure(eq_2d, ic_rare, geom_2d)
            
            # At t = 1.0, the rarefaction fan spans from d = 0 to d = 2.0 along the X axis
            # Y coordinate variations should not affect the planar wave evaluation
            
            @test exact_rare(SVector{3, Float64}(-1.0, 5.0, 1.0))[1] ≈ 0.0 # Left of fan
            @test exact_rare(SVector{3, Float64}(1.0, -3.0, 1.0))[1] ≈ 1.0 # Dead center of the fan (d/t = 1)
            @test exact_rare(SVector{3, Float64}(3.0, 0.0, 1.0))[1]  ≈ 2.0 # Right of fan
        end
    end
end

@testset "Spherical WENO & Alternative Progress Logging" begin
        # 1. Override the global logging behavior to disable the visual ProgressMeter
        # and enforce a tiny interval to guarantee the standard @info logs trigger.
        enable_progress_bar!(false)
        set_progress_interval!(0.001)

        # 2. Build a minimal execution dictionary targeting the untested components
        params = Dict{Symbol, Any}(
            :PDE_name => :linear,
            :PDE_velocities => ((1.0,),),

            :sim_func_name => :run_direct_simulation,
            :tmax => 0.05,
            :snapshots => 2,

            # Hit the Spherical Domain builder
            :Grid_domain => :spherical,
            :Grid_mins => (0.1,),
            :Grid_maxs => (1.0,),
            :Grid_Ns => (250,), # High enough N to guarantee > 1ms execution time for the ETA log
            :Grid_periodic => false,
            :Grid_randomness_factor => (0.0,),
            :Grid_seed => 42,
            :Grid_mover => :none,

            :Weight_name => :exponential,
            :Weight_alpha => 1.0,
            :Weight_range => 2.5,

            :IC_name => :gauss,
            :IC_a => (1.0,),
            :IC_b => (0.5,),
            :IC_width => 0.1,

            :Time_CFL => 0.1,
            :Time_stepper => :RK2,

            # Hit the WENO Scheme builder
            :Scheme_name => :WENO,
            :Scheme_order => 2,
            :Scheme_MLS_order => 0,
            :Flux_name => :Upwind,
            :Limiter_name => :none,

            :MOOD_criterion => :none
        )

        # 3. Execute the simulation and optionally catch the specific logging patterns
        # Using match_mode=:any ignores other logs (like initialization infos)
        @test_logs (:info, r"Using .* for parallel runs!") (:info, r"Simulation Progress: .*") match_mode=:any begin
            sim_data = run_direct_simulation(params)
            
            # Verify the simulation completed its time integration loop normally
            @test sim_data.t[end] ≈ 0.05
            @test length(sim_data.x[end]) == 251
        end
        
        # 4. Restore the default global states so subsequent tests aren't affected
        enable_progress_bar!(true)
        set_progress_interval!(1.0)
    end

using PDEStudio
using GLMakie

@testset "Examples" begin
    @sim advection_2d
    @test true
    @sim burgers_kinetic
    @test true
    @plot euler_1d_sod
    @test true
    reset_plotter!()
end

end