export storage_for, with_backend, set_memory_limit!, get_memory_limit

# =========================================================================
# DEVICE RUNS
#
# The runners place a run where the execution backend of HyperCloudCore runs (`set_backend!`): CPU
# backends keep host arrays, a GPU `KABackend` (e.g. `set_backend!(ROCBackend())` after `using AMDGPU`)
# puts the grid, the scheme and the time stepper on that device via `context[:Storage]`. Setup (initial
# conditions) and output (snapshots) stay on the host and move the first `N` particles in one copy.
# =========================================================================

"""
    storage_for(backend = get_backend()) -> Type

The array type a run lives on for `backend`: `Array` for the CPU backends (`PolyesterBackend`,
`ThreadsBackend`, `KABackend(CPU())`), the device array type of a GPU `KABackend` (e.g. `ROCArray`).
"""
storage_for(::ExecutionBackend = get_backend()) = Array
storage_for(b::KABackend) = Base.typename(typeof(HyperCloudCore.KA.allocate(b.backend, Float32, 0))).wrapper

"""
    with_backend(f, backend)

Calls `f()` with `backend` as the execution backend and restores the previous backend afterwards.
"""
function with_backend(f, backend)
    previous = get_backend()
    set_backend!(backend)
    try
        return f()
    finally
        set_backend!(previous)
    end
end

"""
    build_run_grid(grid_conf, context) -> ParticleGrid

`build_particle_grid` for the storage in `context`. A device grid is constructed on the host first
(neighbour search, volumes) and then moved, so the construction runs on a CPU backend.
"""
function build_run_grid(grid_conf::Dict, context::Dict)
    host_backend = storage_type(context) === Array ? get_backend() : PolyesterBackend()
    return with_backend(() -> build_particle_grid(grid_conf, context), host_backend)
end

"The first `n` entries of a grid array on the host: a view of a host array, a copy of a device array."
host_view(a, n::Integer) = parent(a) isa Array ? view(a, 1:n) : Array(view(a, 1:n))

"""
    set_states!(f, pg::ParticleGrid)

Sets the state of every particle to `f(position)`. Evaluated on the host; a device grid receives the
states in one copy.
"""
function set_states!(f::F, pg::ParticleGrid) where {F}
    N    = pg.meta.N
    pos  = host_view(pg.core.positions, N)
    vals = Vector{eltype(parent(pg.rhos))}(undef, N)
    @inbounds for i in 1:N
        vals[i] = f(pos[i])
    end
    copyto!(pg.rhos, 1, vals, 1, N)
    return nothing
end

# =========================================================================
# MEMORY BUDGET
#
# A machine setting like the backend (not a run parameter, so it does not enter the run hash): with a
# limit, the runners size every growable buffer of the run once, under that many bytes (`reserve!`), and
# a run that outgrows the reservation fails with a `CapacityError` instead of reallocating.
# =========================================================================

const _MEMORY_LIMIT = Ref{Union{Nothing, NamedTuple}}(nothing)

"""
    set_memory_limit!(limit::Integer; tolerance = 0.05, slots_per_particle = nothing)
    set_memory_limit!(nothing)

Reserves the buffers of every following run under `limit` bytes (HyperCloudCore's `reserve!`):
`tolerance` is the share of the limit kept free for the runtime, `slots_per_particle` the neighbour
slots reserved per particle (default: what the initial neighbour lists need; raise it for moving or
adaptive grids). `nothing` removes the limit (buffers grow on demand, the default).
"""
function set_memory_limit!(limit::Integer; tolerance::Real = 0.05, slots_per_particle::Union{Nothing, Real} = nothing)
    limit > 0 || throw(ArgumentError("limit must be positive, got $limit"))
    _MEMORY_LIMIT[] = (; limit = Int(limit), tolerance = Float64(tolerance), slots_per_particle)
end
set_memory_limit!(::Nothing) = (_MEMORY_LIMIT[] = nothing)

"The memory limit set with [`set_memory_limit!`](@ref) as `(; limit, tolerance, slots_per_particle)`, or `nothing`."
get_memory_limit() = _MEMORY_LIMIT[]

"Passes the memory limit to the builders (`context[:MemoryLimit]`, read by `build_timestepper`)."
function apply_memory_limit!(context::Dict)
    m = _MEMORY_LIMIT[]
    m === nothing && return context
    context[:MemoryLimit]     = m.limit
    context[:MemoryTolerance] = m.tolerance
    m.slots_per_particle === nothing || (context[:SlotsPerParticle] = m.slots_per_particle)
    return context
end

"""
    build_run_timestepper(context; time_conf) -> TimeStepper

`build_timestepper`, which also performs the memory reservation. On a device grid the reservation is
verified: with HyperCloudCore 0.6.0 and AMDGPU 2.9.1, `reserve!` can return with zeroed grid buffers when
earlier device runs happened in the same process (allocator reuse, not yet understood). A changed grid
throws instead of silently running on corrupted data.
"""
function build_run_timestepper(context::Dict; time_conf::Dict = Dict{Symbol, Any}())
    checked = haskey(context, :MemoryLimit) && storage_type(context) !== Array
    before  = checked ? _grid_fingerprint(context[:Grid]) : nothing
    ts = build_timestepper(context; time_conf = time_conf)
    if checked && !isequal(_grid_fingerprint(context[:Grid]), before)   # isequal: padding slots may hold NaN
        error("the memory reservation corrupted the device grid (known issue of HyperCloudCore 0.6.0 / " *
              "AMDGPU 2.9.1 after earlier device runs in this process): run in a fresh Julia session or " *
              "without `set_memory_limit!`")
    end
    return ts
end

# Host copies of the first N particles' grid data and of their neighbour slots.
function _grid_fingerprint(pg::ParticleGrid)
    N, c, nb = pg.meta.N, pg.core, pg.neighbor
    ranges = collect(host_view(nb.ranges, N))
    S = maximum(r -> isempty(r) ? 0 : last(r), ranges; init = 0)     # last slot in use (either layout)
    return (ranges, map(a -> collect(host_view(a, N)), (c.positions, c.is_boundary, c.volumes, c.tags, nb.counts)),
            map(a -> collect(host_view(a, S)), (nb.indices, nb.weights, nb.distances)))
end

"Logs the reservation of a run and records it in `sim_data`'s statistics (no-op without a memory limit)."
function record_memory!(sim_data, context::Dict)
    haskey(context, :Memory) || return sim_data
    r = context[:Memory]
    add_stat!(sim_data, :reserved_particles, Float64(r.particles), :none)
    add_stat!(sim_data, :reserved_bytes, Float64(r.bytes), :none)
    return sim_data
end

function _log_memory(context::Dict)
    haskey(context, :Memory) || return
    r = context[:Memory]
    @info "Reserved $(round(r.bytes / 2^20; digits = 1)) MiB of $(round(r.limit / 2^20; digits = 1)) MiB: " *
          "$(r.particles) particles, $(r.slots) neighbour slots"
end
