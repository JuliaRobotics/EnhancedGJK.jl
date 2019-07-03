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
       ReferenceDistance,
       closest_point_in_world,
       closest_point_in_body,
       separation_distance,
       simplex_penetration_distance


include("tagged_points.jl")
include("simplices.jl")
include("johnson_distance.jl")
include("neighbor_mesh.jl")
include("traits.jl")
include("gjk.jl")
include("reference_distance.jl")

end # module
