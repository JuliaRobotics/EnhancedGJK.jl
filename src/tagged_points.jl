
"""
The enhanced GJK algorithm relies on a pre-computed set of neighbors for each
vertex in the mesh. In order to use those neighbors, we have to know from which
vertex to start. Specifically, we need to know the index of the vertex
corresponding to the point in the GJK simplex which we are trying to improve.
To do that, we introduce the notion of a `Tagged` point. A tagged point is
just a point and some arbitrary additional data field. All of the
`any_inside` and `support_vector_max` functions in this package return tagged
points. For most geometries, that tag is empty (`nothing`). But for our
NeighborMesh type, the tag is the linear index into the vertices of the mesh,
which lets us look up that mesh's neighbors faster later on.
"""
struct Tagged{P, T}
    point::P
    tag::T
end

Tagged{P}(point::P) = Tagged(point, nothing)

*(n::Number, t::Tagged) = n * value(t)

value(t::Tagged) = t.point

function any_inside(geometry)
    point = gt.any_inside(geometry)
    Tagged(SVector(point))
end

function any_inside(mesh::gt.AbstractMesh{gt.Point{N, T}}) where {N,T}
    Tagged(SVector(first(gt.vertices(mesh))))
end
