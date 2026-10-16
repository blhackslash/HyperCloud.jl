# Device checks of HyperCloud's physics (PDEs, ICs, builders, kinetic source terms) on HyperCloudCore's
# KABackend(ROCBackend()). Each case mirrors run_direct_simulation / run_kinetic_simulation, but builds the
# grid with context[:Storage] = ROCArray and steps a fixed number of times; the device run must reproduce a
# PolyesterBackend host run. The runners themselves cannot target a device yet (no :Storage plumbing, and
# set_initial_conditions! / saveData! index the grid arrays on the host), so the IC is set on the host here.
using HyperCloud, AMDGPU, StaticArrays, Test, LinearAlgebra
using HyperCloud: get_D, get_M, extract_namespace, build_kinetic_system, set_initial_conditions!, build_ic, build_equation
AMDGPU.allowscalar(false)
const HC = HyperCloud

const DeviceArray    = ROCArray
const device_backend = KABackend(ROCBackend())

# Cases expected to fail, with the reason. Julia reports "Unexpected Pass" once they are fixed.
const KNOWN_BROKEN = Dict(
    r"WENO" => "HyperCloudCore 0.6.0: WENO reconstruction is not kernel-compilable (dynamic basis_length / getindex in WENO.jl:88-94)",
)
is_broken(name) = any(occursin(k, name) for k in keys(KNOWN_BROKEN))

"Build every pipeline component for `params` on `storage`; returns (eq, pg, ts, set_ic!)."
function build_pipeline(params, storage)
    T = get(params, :T, Float64)
    ctx = Dict{Symbol, Any}(:T => T, :Storage => storage)
    ns(k) = extract_namespace(params, k)
    eq_macro = build_equation(ns(:PDE), ctx)
    ctx[:Equation] = eq_macro; ctx[:D] = get_D(eq_macro); ctx[:M] = get_M(eq_macro)
    IC = build_ic(ns(:IC), ctx)
    if params[:sim_func_name] === :run_kinetic_simulation
        eq, st = build_kinetic_system(ns(:Kinetic), ctx)
        ctx[:Equation] = eq; ctx[:M] = get_M(eq)
        ctx[:ExplicitSources] = (); ctx[:ImplicitSources] = (st,)
        set_ic! = pg -> set_initial_conditions!(pg, st, IC, eq_macro)
    else
        eq = eq_macro
        ctx[:ExplicitSources] = (); ctx[:ImplicitSources] = ()
        set_ic! = pg -> set_initial_conditions!(pg, eq, IC)
    end
    ctx[:Domain]     = build_domain(ns(:Grid), ctx)
    ctx[:WeightConf] = ns(:Weight)
    pg = build_particle_grid(ns(:Grid), ctx)
    ctx[:Flux]    = build_flux(ns(:Flux), ctx)
    ctx[:Limiter] = build_limiter(ns(:Limiter), ctx)
    ctx[:MOOD]    = build_mood(ns(:MOOD), ctx)
    ctx[:Scheme]  = build_scheme(ns(:Scheme), ctx)
    ctx[:Tableau] = build_tableau(ns(:Time), ctx)
    ts = build_timestepper(ctx; time_conf = ns(:Time))
    return eq, pg, ts, set_ic!
end

divinterp(ts) = hasproperty(ts, :divergence_interpolator) ? ts.divergence_interpolator : HC.NoFallbackGrad()

"Run `nsteps` on the host and on the device from the same IC and compare."
function compare(params; nsteps = 10)
    T = get(params, :T, Float64)
    set_backend!(PolyesterBackend())                                  # grids are built on the host
    eq, pg_h, ts_h, set_ic! = build_pipeline(params, Array)
    set_ic!(pg_h)
    _, pg_d, ts_d, _ = build_pipeline(params, DeviceArray)
    N = pg_h.meta.N
    @assert pg_d.meta.N == N && Array(pg_d.core.positions)[1:N] == pg_h.core.positions[1:N]
    copyto!(pg_d.rhos, 1, DeviceArray(pg_h.rhos[1:N]), 1, N)
    cfl  = T(get(params, :Time_CFL, 0.2))
    dt_h = cfl * get_time_step(pg_h, eq, divinterp(ts_h))
    set_backend!(device_backend)
    dt_d = cfl * get_time_step(pg_d, eq, divinterp(ts_d))
    for k in 0:nsteps-1; ts_d(eq, pg_d, k * dt_h, dt_h); end
    set_backend!(PolyesterBackend())
    for k in 0:nsteps-1; ts_h(eq, pg_h, k * dt_h, dt_h); end
    ud, uh = Array(pg_d.rhos)[1:N], pg_h.rhos[1:N]
    od, oh = Array(pg_d.core.particle_orders)[1:N], pg_h.core.particle_orders[1:N]
    diff  = maximum(i -> maximum(abs.(ud[i] .- uh[i])), 1:N)
    scale = maximum(i -> maximum(abs.(uh[i])), 1:N)
    return (; diff, scale, dt_h, dt_d, order_mismatches = count(od .!= oh), finite = all(u -> all(isfinite, u), ud))
end

const RESULTS = Pair{String, Any}[]

function check(name, params; nsteps = 10, rtol = 1e-10)
    occursin(Regex(get(ENV, "CASE", "")), name) || return
    @testset "$name" begin
        r = try
            compare(params; nsteps)
        catch e
            push!(RESULTS, name => "ERROR: " * first(sprint(showerror, e), 200))
            is_broken(name) ? (@test_broken false) : (@test false)
            return
        end
        push!(RESULTS, name => r)
        ok = r.finite && isapprox(r.dt_d, r.dt_h; rtol = 1e-12) && r.diff <= rtol * max(r.scale, 1)
        is_broken(name) ? (@test_broken ok) : (@test ok)
        # MOOD orders are not compared strictly: on flat plateaus the detection criterion sees only
        # round-off (~1e-14) and may flip between host and device without changing the state.
    end
end

# --- Configurations --------------------------------------------------------

const KIN_1D = let (v, i) = generate_relax_params(2.0, 1, 1)
    Dict{Symbol, Any}(:Kinetic_velocities => v, :Kinetic_indices => i, :Kinetic_epsilon => 1e-6)
end
const KIN_2D = let (v, i) = generate_relax_params(2.0, 2, 1)
    Dict{Symbol, Any}(:Kinetic_velocities => v, :Kinetic_indices => i, :Kinetic_epsilon => 1e-6)
end

const LIN_1D = merge(Dict{Symbol, Any}(
    :PDE_name => :linear, :PDE_velocities => ((1.0,),), :sim_func_name => :run_direct_simulation,
    :tmax => 0.1, :snapshots => 2,
    :Grid_domain => :rectangular, :Grid_mins => (0.0,), :Grid_maxs => (1.0,), :Grid_Ns => (200,),
    :Grid_periodic => true, :Grid_randomness_factor => (0.0,), :Grid_seed => 42, :Grid_mover => :none,
    :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 5.5,
    :IC_name => :box, :IC_u_bg => (0.1,), :IC_u_box => (1.0,), :IC_mins => (0.3,), :IC_maxs => (0.6,),
    :Time_CFL => 0.2), KIN_1D)

const BURGERS_1D = merge(LIN_1D, Dict{Symbol, Any}(
    :PDE_name => :burgers, :PDE_D => 1, :Grid_mins => (-1.0,), :Grid_maxs => (1.0,),
    :IC_name => :gauss, :IC_a => (1.0,), :IC_b => (0.0,), :IC_width => 0.2))

function lin_2d(periodic::Bool)
    d = merge(Dict{Symbol, Any}(
        :PDE_name => :linear, :PDE_D => 2, :PDE_velocities => ((1.0,), (0.5,)), :sim_func_name => :run_direct_simulation,
        :tmax => 0.1, :snapshots => 2,
        :Grid_domain => :rectangular, :Grid_mins => (0.0, 0.0), :Grid_maxs => (1.0, 1.0), :Grid_Ns => (40, 40),
        :Grid_periodic => periodic, :Grid_randomness_factor => (0.2, 0.2), :Grid_seed => 42, :Grid_mover => :none,
        :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 2.5,
        :IC_name => :box, :IC_u_bg => (0.1,), :IC_u_box => (1.0,), :IC_mins => (0.25, 0.25), :IC_maxs => (0.6, 0.6),
        :Scheme_upwind_alg_nd => :Classic, :MOOD_delta_relax => 1e-6,
        :Time_CFL => 0.2), KIN_2D)
    periodic || (d[:Grid_bc] = Dict{Int, AbstractBoundaryCondition}(i => OutflowBC() for i in 1:4))
    return d
end

const METHODS = sort!([(k, v) for (k, v) in get_master_method_dict() if haskey(v, :Time_stepper)]; by = first)
const SOD     = HyperCloud.euler_1d_sod()

try
    @testset "HyperCloud on $(AMDGPU.device())" begin
        @testset "1D linear advection, box IC" begin
            for (n, m) in METHODS; check("lin1d | $n", merge(LIN_1D, m)); end
        end
        @testset "Sod shock tube (Euler, outflow BC)" begin
            for (n, m) in sort!([(k, v) for (k, v) in SOD.methods if haskey(v, :Time_stepper)]; by = first)
                check("sod | $n", merge(SOD.shared, m))
            end
        end
        @testset "1D Burgers" begin
            check("burgers | RK2 MUSCL2", merge(BURGERS_1D, last(create_method(:RK2, :MUSCL, 2))))
            check("burgers | ARS222 MUSCL2 kinetic", merge(BURGERS_1D, last(create_method(:ARS222, :MUSCL, 2))))
        end
        @testset "2D linear advection (jittered grid)" begin
            for periodic in (true, false), (n, m) in [
                    ("RK2 Upwind1",           last(create_method(:RK2, :Upwind, 1))),
                    ("RK3 MUSCL3 MOOD",       last(create_method(:RK3, :MUSCL, 3; mood_criterion = :U1))),
                    ("ARS222 MUSCL2 kinetic", last(create_method(:ARS222, :MUSCL, 2)))]
                check("2d periodic=$periodic | $n", merge(lin_2d(periodic), m))
            end
        end
        @testset "Float32" begin
            check("f32 | RK2 MUSCL2", merge(LIN_1D, last(create_method(:RK2, :MUSCL, 2)), Dict{Symbol, Any}(:T => Float32)); rtol = 1e-5)
        end
    end
finally
    set_backend!(PolyesterBackend())
    println("\n", "="^30, " results (max |u_dev - u_host|, scale, dt_dev/dt_host - 1, MOOD order mismatches)")
    for (n, r) in RESULTS
        s = r isa String ? r : string("diff=", round(r.diff; sigdigits = 3), "  scale=", round(r.scale; sigdigits = 3),
                                      "  ddt=", round(r.dt_d / r.dt_h - 1; sigdigits = 3), "  mood≠=", r.order_mismatches)
        println(rpad(n, 52), s)
    end
end
