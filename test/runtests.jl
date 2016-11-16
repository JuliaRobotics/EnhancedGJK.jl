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
    result = EnhancedGJK.gjk!(cache, IdentityTransformation(), IdentityTransformation())
    @test isapprox(result.signed_distance, 1.0)
    @test isapprox(result.closest_point_in_body.a, [1.0, 0.0])
    @test isapprox(result.closest_point_in_body.b, [0.0, 0.0])
end

@testset "accelerated mesh gjk" begin
    mesh = DrakeVisualizer.contour_mesh(x -> sum((x - [2, -1.5, 0]).^2) - 1, [-1.5, -1.5, -1.5], [1.5, 1.5, 1.5])
    acc = EnhancedGJK.NeighborMesh(mesh)
    pt = Vec(0., 0, 0)
    cache = EnhancedGJK.CollisionCache(acc, pt)
    result = EnhancedGJK.gjk!(cache, IdentityTransformation(), IdentityTransformation())
    @test isapprox(result.signed_distance, norm([1.2, -0.9, 0]))
    @test isapprox(result.closest_point_in_body.a, [1.2, -0.9, 0.0])
    @test isapprox(result.closest_point_in_body.b, [0.0, 0.0, 0.0])
end

@testset "benchmarks" begin
    include("../perf/runbenchmarks.jl")
end
