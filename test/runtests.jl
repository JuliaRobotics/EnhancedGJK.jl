using EnhancedGJK
import EnhancedGJK: projection_weights, projection_weights_reference
import CoordinateTransformations: IdentityTransformation
import GeometryTypes: Vec
import DrakeVisualizer: contour_mesh
using Base.Test
import StaticArrays: SVector

@testset "johnson distance subalgorithm" begin
    include("johnson_distance.jl")
end

@testset "simplex distance" begin
    simplex = SVector{3}(SVector{2, Float64}[[1., 0], [2., 0], [1., 1]])
    pt = SVector(0., 0)
    cache = EnhancedGJK.CollisionCache(simplex, pt);
    simplex, best_pt, in_interior = EnhancedGJK.gjk!(cache, IdentityTransformation(), IdentityTransformation())
    @show simplex best_pt in_interior
    @test isapprox(norm(best_pt), 1.0)
end

@testset "accelerated mesh gjk" begin
    mesh = DrakeVisualizer.contour_mesh(x -> sum((x - [2, -1.5, 0]).^2) - 1, [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5])
    acc = EnhancedGJK.NeighborMesh(mesh)
    pt = Vec(0., 0, 0)
    cache = EnhancedGJK.CollisionCache(acc, pt)
    simplex, best_pt, in_interior = EnhancedGJK.gjk!(cache, IdentityTransformation(), IdentityTransformation())
    @test isapprox(best_pt, [1.2, -0.9, 0.0])
    @test in_interior == false
end

@testset "benchmarks" begin
    include("../perf/runbenchmarks.jl")
end
