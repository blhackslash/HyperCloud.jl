

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

    @testset "Examples" begin
        @sim advection_1d
        @test true
        @sim burgers_kinetic
        @test true
        @sim euler_1d_sod
        @test true
    end
end
