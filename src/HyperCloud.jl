module HyperCloud

using Reexport

@reexport using PDEStudioCore
@reexport using HyperCloudCore

import HyperCloudCore: flux, max_eigenvalue, prim2cons, cons2prim, velocity

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
include("Domains/_main.jl")
include("SourceTerms/_main.jl")
include("Configs/_main.jl")
include("Builder/_main.jl")

include("inbuilt_methods.jl")
include("run_simulation.jl")
include("time_integration.jl")

end