__precompile__()

module EnhancedGJK

using StaticArrays
import GeometryTypes
const gt = GeometryTypes
using CoordinateTransformations: Transformation,
                                 transform_deriv,
                                 IdentityTransformation
using LinearAlgebra
using Statistics: mean

export CollisionCache,
       gjk!,
       gjk,
       GJKResult,
       NeighborMesh,
       ReferenceDistance


include("tagged_points.jl")
include("simplices.jl")
include("johnson_distance.jl")
include("neighbor_mesh.jl")
include("traits.jl")
include("gjk.jl")
include("reference_distance.jl")

end # module
