module jupiter

using FileIO
using GeometryBasics
using LinearAlgebra
using TensorOperations
using ForwardDiff
using TypedPolynomials

include("meshy.jl")
include("wew.jl")
include("cc.jl")
include("render.jl")

export Mesh
export catmullclark

end