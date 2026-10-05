# =========================================================================
# EXACT CONSTANTS FOR DOUBLE MACH REFLECTION
# =========================================================================
export euler_2d_DMR
# The conservative states: (rho, rho*u, rho*v, E)
# Pre-shock (Right): rho = 1.4, u = 0, v = 0, p = 1.0
const DMR_PRE_SHOCK  = State{4, Float64}(1.4, 0.0, 0.0, 2.5)

# Post-shock (Left): rho = 8.0, u = 7.1447, v = -4.125, p = 116.5
# Energy = 116.5 / 0.4 + 0.5 * 8.0 * (7.1447^2 + (-4.125)^2) = 563.5
const DMR_POST_SHOCK = State{4, Float64}(8.0, 57.15767664977295, -33.0, 563.5)


# =========================================================================
# CUSTOM BOUNDARY CONDITION FUNCTORS
# =========================================================================

# 1. TOP WALL: Time-Dependent Exact Shock Condition
struct DMRTopBC <: AbstractBoundaryCondition end

function (::DMRTopBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    pos = pg.core.positions
    @inbounds for i in 1:pg.meta.N
        if pg.core.tags[i] == tag && pg.core.is_boundary[i]
            x, y = pos[i][1], pos[i][2]
            
            x_shock = 1.0/6.0 + (y + 20.0 * t) / sqrt(3.0)
            rhos_buffer[i] = x < x_shock ? DMR_POST_SHOCK : DMR_PRE_SHOCK
        end
    end
    return nothing
end

# 2. BOTTOM WALL: Slip Wall & Post-Shock Inflow
struct DMRBottomBC <: AbstractBoundaryCondition end

function (::DMRBottomBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    # Step 1: Run the standard OutflowBC to populate the ghost cells
    OutflowBC()(pg, rhos_buffer, tag, ts, eq, t)
    
    pos = pg.core.positions
    @inbounds for i in 1:pg.meta.N
        if pg.core.tags[i] == tag && pg.core.is_boundary[i]
            x = pos[i][1]
            
            if x < 1.0/6.0
                rhos_buffer[i] = DMR_POST_SHOCK
            else
                # Slip wall: flip Y-momentum!
                U_ghost = rhos_buffer[i]
                rhos_buffer[i] = State{4, Float64}(U_ghost[1], U_ghost[2], -U_ghost[3], U_ghost[4])
            end
        end
    end
    return nothing
end

# 3. LEFT WALL: Fixed Inflow
struct DMRLeftBC <: AbstractBoundaryCondition end

function (::DMRLeftBC)(pg, rhos_buffer, tag::Int, ts, eq, t::Real)
    @inbounds for i in 1:pg.meta.N
        if pg.core.tags[i] == tag && pg.core.is_boundary[i]
            rhos_buffer[i] = DMR_POST_SHOCK
        end
    end
    return nothing
end


# =========================================================================
# MAIN SIMULATION SETUP
# =========================================================================

function euler_2d_DMR()
    p0 = (1.0/6.0, 0.0)
    n  = (sqrt(3.0), -1.0) 
    
    bc_map = Dict{Int, AbstractBoundaryCondition}(
        1 => DMRLeftBC(),    # Left Wall
        2 => OutflowBC(),    # Right Wall 
        3 => DMRBottomBC(),  # Bottom Wall
        4 => DMRTopBC(),     # Top Wall
        5 => OutflowBC()     # Corners 
    )

    # 1. Flat Namespace Shared Configuration
    shared = Dict{Symbol, Any}(
        :PDE_name => :euler,
        :PDE_D => 2,
        :PDE_gamma => 1.4,

        :sim_func_name => :run_direct_simulation,                
        :tmax => 0.2,
        :snapshots => 20,
        
        :Grid_domain => :rectangular,
        :Grid_mins => (0.0, 0.0),
        :Grid_maxs => (3.2, 1.0),
        :Grid_Ns => (450, 150),
        :Grid_periodic => false,
        :Grid_bc_map => bc_map,
        :Grid_randomness_factor => (0.2, 0.2),
        :Grid_seed => 42,
        :Grid_mover => :none,
        
        # New Weight Function API mappings
        :Weight_name => :exponential,
        :Weight_alpha => 1.0,
        :Weight_range => 2.5,
        
        :IC_name => :riemann,
        :IC_uL => DMR_POST_SHOCK,
        :IC_uR => DMR_PRE_SHOCK,
        :IC_p0 => p0,
        :IC_n => n,


        :Time_CFL => 0.2,
        
    )

    # 2. Modernized Method Configurations
    methods = Dict{Symbol, Dict{Symbol, Any}}(
        :RK2MUSCL2_EPD2 => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :MUSCL, :Scheme_order => 2, :Scheme_MLS_order => 0, 
            :Flux_name => :Rusanov, :Limiter_name => :none,
            :MOOD_criterion => :U1, :MOOD_strategy => :EPD2, :MOOD_delta_relax => 1e-4
        ),
        :RK2MUSCL2_EPD1 => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :MUSCL, :Scheme_order => 2, :Scheme_MLS_order => 0, 
            :Flux_name => :Rusanov, :Limiter_name => :none,
            :MOOD_criterion => :U1, :MOOD_strategy => :EPD1, :MOOD_delta_relax => 1e-4
        ),
        :RK2MUSCL2_VK => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :MUSCL, :Scheme_order => 2, :Scheme_MLS_order => 0, 
            :Flux_name => :Rusanov, :Limiter_name => :VK, :Limiter_mode => :hard,
            :MOOD_criterion => :none
        ),
        :RK2Upwind1 => Dict(
            :Time_stepper => :RK2, 
            :Scheme_name => :Upwind, :Scheme_order => 1, :Scheme_MLS_order => 1, 
            :Scheme_upwind_alg_nd => :Classic, :Flux_name => :Upwind,
            :MOOD_criterion => :none
        )
    )

    return SimulationConfig(shared, methods, collect(keys(methods)))
end