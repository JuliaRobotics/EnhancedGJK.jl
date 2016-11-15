module EnhancedGJK

using StaticArrays
import StaticArrays: insert
import GeometryTypes
const gt = GeometryTypes
import CoordinateTransformations: Transformation, transform_deriv
import Base: dot, zero, +, *, @pure, convert

@generated function svector{N, T}(v::Union{gt.Vec{N, T}, gt.Point{N, T}})
    svector_impl(v)
end

function svector_impl{N, T}(::Union{Type{gt.Vec{N, T}}, Type{gt.Point{N, T}}})
    Expr(:call, :(SVector{$N, $T}), [:(v[$i]) for i in 1:N]...)
end

include("tagged_points.jl")
include("simplices.jl")
include("johnson_distance.jl")
include("neighbor_mesh.jl")
include("traits.jl")
include("gjk.jl")

end # module
