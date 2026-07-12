module jupiter

using FileIO
using GeometryBasics
using LinearAlgebra
using ForwardDiff
using Serialization
using StaticArrays
using TypedPolynomials

include("meshy.jl")
include("wew.jl")
include("cc.jl")
include("render.jl")
include("reference.jl")

export Mesh
export catmullclark

end