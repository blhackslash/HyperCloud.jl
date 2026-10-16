export storage_for, with_backend

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
