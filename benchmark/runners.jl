# Throughput of the simulation runners: CPU (PolyesterBackend) vs GPU (KABackend(ROCBackend())), Float64 vs
# Float32, on 2-D problems. Uses the GPU test environment:
#     julia -t auto --project=test/gpu benchmark/runners.jl
# BENCH_N=100,200,400,800 (grid points per dimension), BENCH_CASES=<regex>, BENCH_OUT=<file.md>,
# BENCH_GPU=0 (CPU only), BENCH_FTZ=1 (flush subnormals to zero on the CPU workers: Float32 intermediates of
# the neighbour sums can underflow to subnormals, which x86 cores process very slowly and GPUs flush).
# Time per step is the runner's `Runtime` statistic (the time-stepping loop, CFL estimate and snapshots
# included, setup excluded) divided by its step count; every GPU run is compared with the CPU run.
using HyperCloud, AMDGPU, Printf
AMDGPU.allowscalar(false)
enable_progress_bar!(false)

const NS    = parse.(Int, split(get(ENV, "BENCH_N", "100,200,400,800"), ","))
const ONLY  = Regex(get(ENV, "BENCH_CASES", ""))
const OUT   = get(ENV, "BENCH_OUT", "")
const BACKENDS = get(ENV, "BENCH_GPU", "1") == "1" ? [("CPU", PolyesterBackend()), ("GPU", KABackend(ROCBackend()))] :
                                                     [("CPU", PolyesterBackend())]
const FTZ = get(ENV, "BENCH_FTZ", "") == "1"
if FTZ   # every Polyester worker (one per core) and the main thread
    for _ in 1:20
        HyperCloud.HyperCloudCore.@batch per=core for i in 1:10_000; set_zero_subnormals(true); end
    end
    set_zero_subnormals(true)
end

"2-D unit square, `n × n` jittered particles, a run of about 40 CFL steps."
function grid_2d(n; periodic = true)
    d = Dict{Symbol, Any}(
        :tmax => 40 * 0.2 / n, :snapshots => 1, :remove_ghosts => true,
        :Grid_domain => :rectangular, :Grid_mins => (0.0, 0.0), :Grid_maxs => (1.0, 1.0), :Grid_Ns => (n, n),
        :Grid_periodic => periodic, :Grid_randomness_factor => (0.2, 0.2), :Grid_seed => 42, :Grid_mover => :none,
        :Weight_name => :exponential, :Weight_alpha => 1.0, :Weight_range => 2.5,
        :Time_CFL => 0.2, :Scheme_upwind_alg_nd => :Classic, :MOOD_delta_relax => 1e-4)
    periodic || (d[:Grid_bc] = Dict{Int, AbstractBoundaryCondition}(i => OutflowBC() for i in 1:4))
    return d
end

const ADVECTION = Dict{Symbol, Any}(:PDE_name => :linear, :PDE_D => 2, :PDE_velocities => ((1.0,), (0.5,)),
    :IC_name => :box, :IC_u_bg => (0.1,), :IC_u_box => (1.0,), :IC_mins => (0.25, 0.25), :IC_maxs => (0.6, 0.6))
const EULER = Dict{Symbol, Any}(:PDE_name => :euler, :PDE_D => 2, :PDE_gamma => 1.4,
    :IC_name => :riemann, :IC_uL => (1.0, 0.0, 0.0, 2.5), :IC_uR => (0.125, 0.0, 0.0, 0.25),
    :IC_p0 => (0.5, 0.5), :IC_n => (1.0, 1.0))                     # diagonal Sod tube (conservative states)
const KINETIC = let (v, i) = generate_relax_params(2.0, 2, 1)
    Dict{Symbol, Any}(:Kinetic_velocities => v, :Kinetic_indices => i, :Kinetic_epsilon => 1e-6)
end

const CASES = [
    ("A advection RK2 MUSCL2",        n -> merge(grid_2d(n), ADVECTION, last(create_method(:RK2, :MUSCL, 2)))),
    ("B advection RK3 MUSCL3 MOOD",   n -> merge(grid_2d(n), ADVECTION, last(create_method(:RK3, :MUSCL, 3; mood_criterion = :U1)))),
    ("C Euler RK2 MUSCL2 MOOD",       n -> merge(grid_2d(n; periodic = false), EULER, last(create_method(:RK2, :MUSCL, 2; mood_criterion = :U1)))),
    ("D kinetic ARS222 MUSCL2",       n -> merge(grid_2d(n), ADVECTION, KINETIC, last(create_method(:ARS222, :MUSCL, 2)))),
]

stat(sim, k) = sim.stats[k][1]

function run_on(backend, params)
    run = getfield(HyperCloud, params[:sim_func_name])
    sim = with_backend(() -> run(params), backend)
    sim === nothing && error("run failed: $(params[:sim_func_name]) on $backend")
    return sim
end

maxdiff(a, b) = maximum(i -> maximum(abs.(a[i] .- b[i])), eachindex(a))
scale(a)      = maximum(u -> maximum(abs, u), a)

rows = NamedTuple[]
for (name, make) in CASES
    occursin(ONLY, name) || continue
    for T in (Float64, Float32)
        for (_, b) in BACKENDS; run_on(b, merge(make(20), Dict{Symbol, Any}(:T => T))); end   # compile
        for n in NS
            params = merge(make(n), Dict{Symbol, Any}(:T => T))
            sims = Dict(bn => run_on(b, params) for (bn, b) in BACKENDS)
            for (bn, _) in BACKENDS
                s = sims[bn]
                steps, rt, N = Int(stat(s, :k_step)), stat(s, :Runtime), Int(stat(s, :final_particle_count))
                push!(rows, (; case = name, T = string(T), n, N, backend = bn, steps, ms_step = 1e3 * rt / steps,
                               mpps = N * steps / rt / 1e6,
                               rel_diff = bn == "GPU" ? maxdiff(sims["GPU"].u[end], sims["CPU"].u[end]) / scale(sims["CPU"].u[end]) : NaN,
                               same_steps = bn == "GPU" ? stat(s, :k_step) == stat(sims["CPU"], :k_step) : true))
                r = rows[end]
                @printf("%-30s %-7s %4d² %-3s %4d steps  %9.3f ms/step  %8.2f Mpart-steps/s  %s\n", name, r.T, n, bn,
                        steps, r.ms_step, r.mpps, bn == "GPU" ? @sprintf("rel diff to CPU %.1e%s", r.rel_diff, r.same_steps ? "" : "  STEP COUNT DIFFERS") : "")
            end
        end
    end
end

# Markdown summary: one row per case, precision and size; speed-up = CPU time / GPU time per step.
io = IOBuffer()
if length(BACKENDS) == 1
    println(io, "| case | T | particles | CPU ms/step | CPU Mpart-steps/s |\n|---|---|---:|---:|---:|")
    foreach(r -> @printf(io, "| %s | %s | %d | %.2f | %.2f |\n", r.case, r.T, r.N, r.ms_step, r.mpps), rows)
else
    println(io, "| case | T | particles | CPU ms/step | GPU ms/step | GPU speed-up | GPU rel. diff |")
    println(io, "|---|---|---:|---:|---:|---:|---:|")
end
for r in rows
    r.backend == "GPU" || continue
    c = only(filter(x -> x.case == r.case && x.T == r.T && x.n == r.n && x.backend == "CPU", rows))
    @printf(io, "| %s | %s | %d | %.2f | %.2f | %.1f× | %.0e |\n", r.case, r.T, r.N, c.ms_step, r.ms_step, c.ms_step / r.ms_step, r.rel_diff)
end
table = String(take!(io))
println("\n", table)
isempty(OUT) || write(OUT, "CPU: $(Sys.cpu_info()[1].model), $(Threads.nthreads()) Julia threads (Polyester: one worker per core)" *
                         (FTZ ? ", subnormals flushed to zero" : "") * "\n" *
                         "GPU: $(AMDGPU.device())\n\n" * table)
