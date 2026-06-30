module jupiter

using FileIO
using GeometryBasics
using LinearAlgebra
using TensorOperations
using TaylorDiff
using TypedPolynomials

include("meshy.jl")
include("wew.jl")
include("cc.jl")

export Mesh
export catmullclark

end