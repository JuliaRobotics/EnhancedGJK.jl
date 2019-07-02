using EnhancedGJK
using MeshCat
using Polyhedra
using StaticArrays: SVector
using GeometryTypes: HyperSphere, Point

function visualize_simplex(vis::Visualizer, simplex)
    p = polyhedron(vrep(simplex))
    setobject!(vis[:simplex], Polyhedra.Mesh(p))
    for (i, point) in enumerate(simplex)
        setobject!(vis["p$i"], HyperSphere(Point(point), 0.03))
    end
end

vis = Visualizer()
open(vis)
simplex = SVector{4, SVector{3, Float64}}([-3.33728783796428, 0.3321305518800686, 0.11004228580261734], [-3.33728783796428, -0.6678694481199314, 0.11004228580261734], [1.6627121620357201, -0.6678694481199314, 0.11004228580261734], [1.6627121620357201, 0.3321305518800686, 0.11004228580261734])
visualize_simplex(vis, simplex)
@show EnhancedGJK.projection_weights(simplex)
