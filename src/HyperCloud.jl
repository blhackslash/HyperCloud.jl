module HyperCloud

using Reexport

@reexport using PDEStudioCore
@reexport using HyperCloudCore

import HyperCloudCore: flux, max_eigenvalue, prim2cons, cons2prim, velocity, implicit_solve, evaluate_source, pre_solve_update!

using StaticArrays
using LinearAlgebra
using Random
using ProgressMeter

const _PROGRESS_BAR = Ref(true)
const _PROGRESS_INTERVAL= Ref(1.)

enable_progress_bar!(b::Bool) = (_PROGRESS_BAR[] = b)
set_progress_interval!(t::Real) = (_PROGRESS_INTERVAL[] = Float64(t))

include("ICs/_main.jl")
include("PDEs/_main.jl")
include("SourceTerms/_main.jl")
include("Examples/_main.jl")
include("Builder/_main.jl")
include("SimulationFunctions/_main.jl")

include("inbuilt_methods.jl")
include("time_integration.jl")
include("macros.jl")

export @sim, @plot, @backup
export enable_progress_bar!, set_progress_interval!
export BONITO_LISTEN_URL, BONITO_LISTEN_PORT, BONITO_PROXY_URL, BONITO_IS_CONFIGURED
export CONFIG_SEARCH_PATHS, CONFIG_ROOT_OVERWRITE

function __init__()
    set_target_module!(@__MODULE__)
end
end