module HyperCloud

using StaticArrays
using PDEStudioCore
using HyperCloudCore
using LinearAlgebra

include("ICs/_main.jl")
include("PDEs/_main.jl")
include("Domains/_main.jl")
include("SourceTerms/_main.jl")

include("inbuilt_methods.jl")

end