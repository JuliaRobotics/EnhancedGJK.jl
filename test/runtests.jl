using Test
using Statistics
using LinearAlgebra
using Random

using EnhancedGJK
using EnhancedGJK: projection_weights, projection_weights_reference
using CoordinateTransformations: IdentityTransformation, Translation
using StaticArrays: SVector
import GeometryTypes
const gt = GeometryTypes
using FileIO

const mesh_dir = joinpath(dirname(@__FILE__), "meshes")

@testset "reference distance" begin
    mesh = load(joinpath(mesh_dir, "base_link.obj"))
    for x in range(0.05, stop=1, length=10)
        for y in range(-1, stop=1, length=10)
            for z in range(-1, stop=1, length=10)
                point = SVector(x, y, z)
                gjk_dist = gjk(mesh, point)
                simple_dist = ReferenceDistance.signed_distance(mesh, point)
                @test isapprox(gjk_dist.signed_distance, simple_dist, atol=2e-4)
            end
        end
    end
end

@testset "reference interior distance" begin
    verts = [gt.Point(1., -1, 1), gt.Point(1., 1, 1), gt.Point(-1., 1, 1), gt.Point(-1., -1, 1),
             gt.Point(1., -1, -1), gt.Point(1., 1, -1), gt.Point(-1., 1, -1), gt.Point(-1., -1, -1)]
    ft = gt.Face{3, Int}
    faces = [ft(1, 2, 5), ft(6, 5, 2),
             ft(2, 3, 6), ft(7, 6, 3),
             ft(3, 4, 7), ft(8, 7, 4),
             ft(4, 1, 8), ft(5, 8, 1),
             ft(2, 1, 3), ft(4, 3, 1),
             ft(5, 6, 7), ft(7, 8, 5)]
    mesh = gt.HomogenousMesh(verts, faces)
    ReferenceDistance.signed_distance(mesh, SVector(0, 0, 0))

    for x in range(-1, stop=1, length=10)
        for y in range(-1, stop=1, length=10)
            for z in range(-1, stop=1, length=10)
                point = SVector(x, y, z)
                expected = -min(1 - abs(x), 1 - abs(y), 1 - abs(z))
                actual = ReferenceDistance.signed_distance(mesh, point)
                @test isapprox(expected, actual, atol=1e-15)
            end
        end
    end
end

@testset "table" begin
    width = 0.5
    thickness = 0.05
    surface_points = Vector{SVector{3, Float64}}()
    for z in [-thickness, thickness]
    for x in [-width, width]
        for y in [-width, width]
            push!(surface_points, SVector(x, y, z))
        end
    end
    end
    geometry = SVector{length(surface_points)}(surface_points)
    point = zeros(SVector{3, Float64})
    for x in range(-width, stop=width, length=21)
        for y in range(-width, stop=width, length=21)
            z = 0.1
            @test isapprox(gjk(geometry,
                point,
                IdentityTransformation(),
                Translation(SVector(x, y, z))).signed_distance, 0.05)
            z = 0.06
            @test isapprox(gjk(geometry,
                    point,
                    IdentityTransformation(),
                    Translation(SVector(x, y, z))).signed_distance,
                0.01,
                atol=1e-12)
            z = 0.05
            @test isapprox(gjk(geometry,
                    point,
                    IdentityTransformation(),
                    Translation(SVector(x, y, z))).signed_distance,
                0.0,
                atol=1e-12)
        end
    end
end



@testset "johnson distance subalgorithm" begin
    include("johnson_distance.jl")
end

@testset "simplex distance" begin
    simplex = SVector{3}(SVector{2, Float64}[[1., 0], [2., 0], [1., 1]])
    pt = SVector(0., 0)
    cache = CollisionCache(simplex, pt);
    result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
    @test isapprox(result.signed_distance, 1.0)
    @test isapprox(result.closest_point_in_body.a, [1.0, 0.0])
    @test isapprox(result.closest_point_in_body.b, [0.0, 0.0])
end

@testset "mesh to mesh" begin
    mesh = load(joinpath(mesh_dir, "r_foot_chull.obj"))
    dx = 1.0
    foot_length = 0.172786 + 0.090933

    cache = CollisionCache(mesh, mesh)
    result = gjk!(cache, IdentityTransformation(), Translation(SVector(dx, 0, 0)))
    @test isapprox(result.signed_distance, dx - foot_length, atol=1e-3)

    cache = CollisionCache(mesh, mesh)
    result = gjk!(cache, Translation(SVector(dx, 0, 0)), IdentityTransformation())
    @test isapprox(result.signed_distance, dx - foot_length, atol=1e-3)

    cache = CollisionCache(mesh, mesh)
    expected_penetration = 0.01
    result = gjk!(cache, IdentityTransformation(), Translation(foot_length - expected_penetration, 0, 0))
    # TODO: penetration distance is inconsistent and inaccurate
    @test result.signed_distance < 0
end

@testset "neighbor mesh to mesh" begin
    mesh = NeighborMesh(load(joinpath(mesh_dir, "r_foot_chull.obj")))
    dx = 1.0
    foot_length = 0.172786 + 0.090933

    cache = CollisionCache(mesh, mesh)
    result = gjk!(cache, IdentityTransformation(), Translation(SVector(dx, 0, 0)))
    @test isapprox(result.signed_distance, dx - foot_length, atol=1e-3)

    cache = CollisionCache(mesh, mesh)
    result = gjk!(cache, Translation(SVector(dx, 0, 0)), IdentityTransformation())
    @test isapprox(result.signed_distance, dx - foot_length, atol=1e-3)

    cache = CollisionCache(mesh, mesh)
    expected_penetration = 0.01
    result = gjk!(cache, IdentityTransformation(), Translation(foot_length - expected_penetration, 0, 0))
    # TODO: penetration distance is inconsistent and inaccurate
    @test result.signed_distance < 0
end

@testset "geometry types" begin
    # Adapted from
    # https://github.com/JuliaGeometry/GeometryTypes.jl/blob/master/test/gjk.jl
    @testset "gjk examples" begin
        c1 = gt.Simplex(gt.Vec(-1.))
        c2 = gt.Simplex(gt.Vec(4.))
        cache = CollisionCache(c1, c2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(result.signed_distance, 5.0)

        c1 = gt.Simplex(gt.Vec(-1.,0,0))
        c2 = gt.Simplex(gt.Vec(4.,0,0))
        cache = CollisionCache(c1, c2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(result.signed_distance, 5.0)

        c1 = gt.FlexibleConvexHull([gt.Vec(0.,0), gt.Vec(0.,1), gt.Vec(1.,0),gt.Vec(1.,1)])
        c2 = gt.Simplex(gt.Vec(4.,0.5))
        cache = CollisionCache(c1, c2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(result.signed_distance, 3.0)

        pt1 = gt.Vec(1,2,3.)
        pt2 = gt.Vec(3,4,5.)
        cache = CollisionCache(pt1, pt2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(result.signed_distance, norm(pt1 - pt2))
    end

    @testset "gjk intersecting lines" begin
        c1 = gt.Simplex(gt.Vec(1,1.), gt.Vec(1, 2.))
        @test gjk(c1, c1).signed_distance == 0.

        c2 = gt.Simplex(gt.Vec(1,1.), gt.Vec(10, 2.))
        @test gjk(c1, c2).signed_distance == 0.
    end
end

@testset "benchmarks" begin
    include("../perf/runbenchmarks.jl")
end
