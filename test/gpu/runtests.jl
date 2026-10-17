# Device checks of HyperCloud on HyperCloudCore's KABackend(ROCBackend()); every device run must reproduce a
# PolyesterBackend host run.
#  - "pipeline": the runners' build steps with context[:Storage] = ROCArray, a fixed number of steps
#    (isolates the physics: PDEs, ICs, builders, kinetic source terms);
#  - "runner": complete run_direct_simulation / run_kinetic_simulation runs, snapshots compared.
using HyperCloud, AMDGPU, StaticArrays, Test, LinearAlgebra, Logging
using HyperCloud: get_D, get_M, extract_namespace, build_kinetic_system, set_initial_conditions!, build_ic, build_equation
AMDGPU.allowscalar(false)
const HC = HyperCloud

const DeviceArray    = ROCArray
const device_backend = KABackend(ROCBackend())

# Cases expected to fail, with the reason. Julia reports "Unexpected Pass" once they are fixed. (Empty since the
# WENO kernel and the :sinusoidal mover were fixed in HyperCloudCore after 0.6.0.)
const KNOWN_BROKEN = Dict{Regex, String}()
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
    pg = HC.build_run_grid(ns(:Grid), ctx)
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
    _, pg_d, ts_d, set_ic_d! = build_pipeline(params, DeviceArray)
    set_ic_d!(pg_d)
    N = pg_h.meta.N
    @assert pg_d.meta.N == N && Array(pg_d.core.positions)[1:N] == pg_h.core.positions[1:N]
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

"Run `params` with the runner on the host and on the device and compare the snapshots."
function compare_runner(params)
    run = getfield(HyperCloud, params[:sim_func_name])
    h = with_backend(() -> run(params), PolyesterBackend())
    d = with_backend(() -> run(params), device_backend)
    (h === nothing || d === nothing) && error("run failed (host: $(h !== nothing), device: $(d !== nothing))")
    length(h.t) == length(d.t) || error("snapshot counts differ: $(length(h.t)) vs $(length(d.t))")
    maxdiff(a, b) = maximum(k -> isempty(a[k]) ? 0.0 : maximum(i -> maximum(abs.(a[k][i] .- b[k][i])), eachindex(a[k])), eachindex(a))
    all(k -> length(h.x[k]) == length(d.x[k]), eachindex(h.x)) || error("particle counts differ")
    return (; diff = maxdiff(h.u, d.u), scale = maximum(k -> maximum(u -> maximum(abs, u), h.u[k]), eachindex(h.u)),
              xdiff = maxdiff(h.x, d.x), dt_h = h.t[end], dt_d = d.t[end], nsnap = length(h.t),
              N = length(h.x[end]), order_mismatches = 0, finite = all(k -> all(u -> all(isfinite, u), d.u[k]), eachindex(d.u)))
end

function check_runner(name, params; rtol = 1e-10)
    occursin(Regex(get(ENV, "CASE", "")), name) || return
    @testset "$name" begin
        r = try
            compare_runner(params)
        catch e
            push!(RESULTS, name => "ERROR: " * first(sprint(showerror, e), 200))
            is_broken(name) ? (@test_broken false) : (@test false)
            return
        end
        push!(RESULTS, name => r)
        ok = r.finite && r.dt_d == r.dt_h && r.diff <= rtol * max(r.scale, 1) && r.xdiff <= rtol
        is_broken(name) ? (@test_broken ok) : (@test ok)
    end
end

const METHODS = sort!([(k, v) for (k, v) in get_master_method_dict() if haskey(v, :Time_stepper)]; by = first)
const SOD     = HyperCloud.euler_1d_sod()

try
    @testset "HyperCloud on $(AMDGPU.device())" begin
        # First, before any other device run: with HyperCloudCore 0.6.0 / AMDGPU 2.9.1 a reservation after
        # earlier device runs in the process can zero grid buffers (the runners detect it and refuse to run).
        @testset "Memory limit" begin
            enable_progress_bar!(false)
            try
                set_memory_limit!(2^30; slots_per_particle = 19)            # reserved once, nothing grows
                sod_adapt = merge(SOD.shared, SOD.methods[:RK2MUSCL2_U1_EPD1], Dict{Symbol, Any}(:Grid_mover => :lagrangian, :Adapt_enabled => true))
                check_runner("runner | sod lagrangian mover + adaptivity, memory limit", sod_adapt)
                # A further reserved device run in the same process: correct, or refused by the runner's check.
                @testset "second reserved device run: correct or refused" begin
                    logs, sim = Test.collect_test_logs(min_level = Logging.Error) do
                        with_backend(() -> run_direct_simulation(merge(lin_2d(false), last(create_method(:ARS222, :MUSCL, 2)))), device_backend)
                    end
                    refused = any(l -> occursin("memory reservation corrupted", sprint(showerror, first(l.kwargs[:exception]))), logs)
                    push!(RESULTS, "runner | second reserved device run" => (sim !== nothing ? "ran" : refused ? "refused by the corruption check" : "FAILED otherwise"))
                    @test sim !== nothing || refused
                end
            finally
                set_memory_limit!(nothing)
            end
        end
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
        @testset "Runners" begin
            enable_progress_bar!(false)
            mood = Dict{Symbol, Any}(:MOOD_delta_relax => 1e-6)
            check_runner("runner | lin1d RK3 MUSCL3 MOOD", merge(LIN_1D, last(create_method(:RK3, :MUSCL, 3; mood_criterion = :U1)), mood))
            check_runner("runner | lin1d ARS222 MUSCL2 kinetic", merge(LIN_1D, last(create_method(:ARS222, :MUSCL, 2))))
            check_runner("runner | lin1d ARS222 MUSCL2 kinetic save_relax", merge(LIN_1D, last(create_method(:ARS222, :MUSCL, 2)), Dict{Symbol, Any}(:save_relax => true)))
            for (n, m) in sort!([(k, v) for (k, v) in SOD.methods if haskey(v, :Time_stepper)]; by = first)
                check_runner("runner | sod $n", merge(SOD.shared, m))
            end
            check_runner("runner | sod RK2Upwind1 keep ghosts", merge(SOD.shared, SOD.methods[:RK2Upwind1], Dict{Symbol, Any}(:remove_ghosts => false)))
            for periodic in (true, false)
                check_runner("runner | 2d periodic=$periodic RK3 MUSCL3 MOOD", merge(lin_2d(periodic), last(create_method(:RK3, :MUSCL, 3; mood_criterion = :U1))))
                check_runner("runner | 2d periodic=$periodic ARS222 MUSCL2 kinetic", merge(lin_2d(periodic), last(create_method(:ARS222, :MUSCL, 2))))
            end
            lagr = Dict{Symbol, Any}(:Grid_mover => :lagrangian)
            check_runner("runner | burgers lagrangian mover", merge(BURGERS_1D, last(create_method(:RK2, :MUSCL, 2)), lagr))
            check_runner("runner | sod lagrangian mover", merge(SOD.shared, SOD.methods[:RK2MUSCL2_U1_EPD1], lagr))
            check_runner("runner | sod lagrangian mover + adaptivity", merge(SOD.shared, SOD.methods[:RK2MUSCL2_U1_EPD1], lagr, Dict{Symbol, Any}(:Adapt_enabled => true)))
            check_runner("runner | sod scaled lagrangian mover + adaptivity", merge(SOD.shared, SOD.methods[:RK2MUSCL2_U1_EPD1],
                         Dict{Symbol, Any}(:Grid_mover => :scaled_lagrangian, :Grid_mover_alpha => 0.5, :Adapt_enabled => true)))
            check_runner("runner | lin1d sinusoidal mover", merge(LIN_1D, last(create_method(:RK2, :MUSCL, 2)), Dict{Symbol, Any}(:Grid_mover => :sinusoidal, :Grid_mover_amplitude => 0.02)))
            check_runner("runner | 2d MLS volumes", merge(lin_2d(true), last(create_method(:RK2, :MUSCL, 2)), Dict{Symbol, Any}(:Grid_volumes => :mls)))
            check_runner("runner | 2d MLS volumes, sinusoidal mover", merge(lin_2d(true), last(create_method(:RK2, :MUSCL, 2)),
                         Dict{Symbol, Any}(:Grid_volumes => :mls, :Grid_mover => :sinusoidal, :Grid_mover_amplitude => 0.02)))
            check_runner("runner | f32 lin1d RK2 MUSCL2", merge(LIN_1D, last(create_method(:RK2, :MUSCL, 2)), Dict{Symbol, Any}(:T => Float32)); rtol = 1e-5)
            check_runner("runner | lin1d RK2WENO2", merge(LIN_1D, get_master_method_dict()["RK2WENO2"]))
            enable_progress_bar!(true)
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
                                      "  ddt=", round(r.dt_d / r.dt_h - 1; sigdigits = 3),
                                      haskey(r, :nsnap) ? "  xdiff=$(round(r.xdiff; sigdigits = 3))  N=$(r.N)" : "  mood≠=$(r.order_mismatches)")
        println(rpad(n, 52), s)
    end
end
