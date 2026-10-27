export bgk_1d_sod

"""
    bgk_1d_sod()

The Sod shock tube for the discrete-velocity BGK model (`γ = 1.4`, i.e. a Chu-reduced distribution
pair on a 1-D velocity grid of 24 nodes, 48 equations) at the relaxation times `τ = 1e-5, 1e-3, 1e-2`,
against the exact Euler solution ([`euler_limit_solution`](@ref)), which the BGK solution approaches
as `τ → 0`. Methods: IMEX `ARS222` with MUSCL2 (Venkatakrishnan limiter, or U2 MOOD) and explicit
`RK2` with MUSCL2 (the collision term explicit, only stable for the larger `τ`).
"""
function bgk_1d_sod()
    varied_params = create_varied_dict()
    push!(varied_params, :PDE_tau => [1e-5, 1e-3, 1e-2])

    bc_map = Dict{Int, AbstractBoundaryCondition}(
        1 => FixedDirichlet(), # Left Boundary
        2 => OutflowBC()       # Right Boundary
    )

    shared = Dict{Symbol, Any}(
        # Discrete-velocity BGK model; the IC below is macroscopic (Euler conservative variables)
        :PDE_name => :bgk,
        :PDE_D => 1,
        :PDE_gamma => 1.4,
        :PDE_Nv => 24,
        :PDE_vmax => 6.0,
        :PDE_velocity_grid => :uniform,
        :PDE_tau => 1e-5,

        :sim_func_name => :run_bgk_simulation,
        :tmax => 0.2,
        :snapshots => 20,

        :Grid_domain => :rectangular,
        :Grid_mins => (-0.5,),
        :Grid_maxs => (0.5,),
        :Grid_Ns => (330,),
        :Grid_periodic => false,
        :Grid_bc => bc_map,
        :Grid_randomness_factor => (0.2,),
        :Grid_seed => 10,
        :Grid_mover => :none,

        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 5.5,

        :IC_name => :riemann,
        :IC_uL => SOD_LEFT,
        :IC_uR => SOD_RIGHT,
        :IC_p0 => (0.0,),
        :IC_n => (1.0,),

        :Time_CFL => 0.2,
    )

    methods = Dict{Symbol, Dict{Symbol, Any}}(
        :exact_euler_limit => Dict(),
        :ARS222MUSCL2_VK => Dict(
            :Time_stepper => :ARS222,
            :Scheme_name => :MUSCL, :Scheme_order => 2,
            :Flux_name => :Rusanov, :Limiter_name => :VK, :Limiter_mode => :hard,
            :MOOD_criterion => :none
        ),
        :ARS222MUSCL2_U2_EPD1 => Dict(
            :Time_stepper => :ARS222,
            :Scheme_name => :MUSCL, :Scheme_order => 2,
            :Flux_name => :Rusanov, :Limiter_name => :none,
            :MOOD_criterion => :U2, :MOOD_strategy => :EPD1, :MOOD_delta_relax => 1e-4
        ),
        :RK2MUSCL2_VK => Dict(
            :Time_stepper => :RK2,
            :Scheme_name => :MUSCL, :Scheme_order => 2,
            :Flux_name => :Rusanov, :Limiter_name => :VK, :Limiter_mode => :hard,
            :MOOD_criterion => :none
        ),
    )

    active_methods = [:exact_euler_limit, :ARS222MUSCL2_VK, :ARS222MUSCL2_U2_EPD1]

    return SimulationConfig(shared, methods, active_methods; varied = varied_params, reference = :euler_limit_solution)
end
