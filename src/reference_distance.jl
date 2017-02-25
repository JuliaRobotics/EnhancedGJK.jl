module ReferenceDistance

using GeometryTypes
import EnhancedGJK: svector, projection_weights
import Base: convert
import StaticArrays: SVector

immutable Plane{N, T}
    a::SVector{N, T}
    b::T
end

function normal(points)
    span = GeometryTypes.edgespan(Simplex(points))
    normalize(cross(Point(column(span, 1)), Point(column(span, 2))))
end

function planefit{N, T}(points::NTuple{3, Point{N, T}})
    a = normal(points)
    b = -dot(a, points[1])
    Plane(svector(a), b)
end

@generated function convert{S <: SVector, N}(::Type{S}, simplex::Simplex{N})
    Expr(:call, :SVector, Expr(:tuple, [:(svector(simplex[$i])) for i in 1:N]...))
end

function interior_distance(face_points, target)
    plane = planefit(face_points)
    dot(plane.a, target) + plane.b
end

function exterior_distance(face_points, target)
    simplex = convert(SVector, Simplex(face_points)) .- SVector((target,))
    weights = projection_weights(simplex)
    projected = dot(weights, simplex)
    norm(projected)
end

# function reverse{T, I}(face::Face{3, T, I})
#     Face{3, T, I}(face[3], face[2], face[2])
# end

# function oriented_faces(mesh)
#     fs = faces(mesh)
#     vs = vertices(mesh)
#     ns = normals(mesh)
#     for (i, f) in enumerate(fs)
#         avg_normal = mean(ns[f])
#         if dot(avg_normal, normal(vs[f])) < 0
#             fs[i] = reverse(fs[i])
#         end
#     end
#     fs
# end

function signed_distance(mesh::AbstractMesh, point)
    verts = vertices(mesh)
    fs = faces(mesh)
    int_dist = minimum([interior_distance(verts[f], point) for f in fs])
    if int_dist >= 0
        return -int_dist
    else
        ext_dist = minimum([exterior_distance(verts[f], point) for f in fs])
        return ext_dist
    end
end

function signed_distance(mesh::AbstractMesh)
    point -> signed_distance(mesh, point)
end

end
