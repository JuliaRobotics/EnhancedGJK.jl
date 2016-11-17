__precompile__()

module EnhancedGJK

using StaticArrays
import GeometryTypes
const gt = GeometryTypes
import CoordinateTransformations: Transformation,
                                  transform_deriv,
                                  IdentityTransformation
import Base: dot,
             zero,
             +,
             *,
             @pure,
             convert

export CollisionCache,
       gjk!,
       gjk,
       GJKResult,
       NeighborMesh


include("geometry_type_conversions.jl")
include("tagged_points.jl")
include("simplices.jl")
include("johnson_distance.jl")
include("neighbor_mesh.jl")
include("traits.jl")
include("gjk.jl")

end # module
