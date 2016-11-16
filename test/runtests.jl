using EnhancedGJK
import EnhancedGJK: projection_weights, projection_weights_reference
import CoordinateTransformations: IdentityTransformation, Translation
import GeometryTypes: Vec
using FileIO
using MeshIO
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

@testset "mesh to mesh" begin
    mesh = load("meshes/r_foot_chull.obj")
    dx = 1.0
    foot_length = 0.172786 + 0.090933

    cache = EnhancedGJK.CollisionCache(mesh, mesh)
    result = EnhancedGJK.gjk!(cache, IdentityTransformation(), Translation(SVector(dx, 0, 0)))
    @test isapprox(result.signed_distance, dx - foot_length, atol=1e-3)

    cache = EnhancedGJK.CollisionCache(mesh, mesh)
    result = EnhancedGJK.gjk!(cache, Translation(SVector(dx, 0, 0)), IdentityTransformation())
    @test isapprox(result.signed_distance, dx - foot_length, atol=1e-3)

    cache = EnhancedGJK.CollisionCache(mesh, mesh)
    expected_penetration = 0.01
    result = EnhancedGJK.gjk!(cache, IdentityTransformation(), Translation(foot_length - expected_penetration, 0, 0))
    # TODO: penetration distance is inconsistent and inaccurate
    @test result.signed_distance < 0
end

@testset "neighbor mesh to mesh" begin
    mesh = EnhancedGJK.NeighborMesh(load("meshes/r_foot_chull.obj"))
    dx = 1.0
    foot_length = 0.172786 + 0.090933

    cache = EnhancedGJK.CollisionCache(mesh, mesh)
    result = EnhancedGJK.gjk!(cache, IdentityTransformation(), Translation(SVector(dx, 0, 0)))
    @test isapprox(result.signed_distance, dx - foot_length, atol=1e-3)

    cache = EnhancedGJK.CollisionCache(mesh, mesh)
    result = EnhancedGJK.gjk!(cache, Translation(SVector(dx, 0, 0)), IdentityTransformation())
    @test isapprox(result.signed_distance, dx - foot_length, atol=1e-3)

    cache = EnhancedGJK.CollisionCache(mesh, mesh)
    expected_penetration = 0.01
    result = EnhancedGJK.gjk!(cache, IdentityTransformation(), Translation(foot_length - expected_penetration, 0, 0))
    # TODO: penetration distance is inconsistent and inaccurate
    @test result.signed_distance < 0
end

@testset "benchmarks" begin
    include("../perf/runbenchmarks.jl")
end
