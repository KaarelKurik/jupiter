module jupiter

using FileIO
using GeometryBasics
using LinearAlgebra
using ForwardDiff
using Serialization
using StaticArrays
using TypedPolynomials

include("meshy.jl")
include("cc.jl")
include("throat.jl")
include("ad.jl")
include("chart.jl")
include("geodesic.jl")
include("mouth.jl")
include("fit.jl")
include("render.jl")
include("reference.jl")

export Mesh
export catmullclark

end