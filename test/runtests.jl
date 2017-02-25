using EnhancedGJK
import EnhancedGJK: projection_weights, projection_weights_reference
import CoordinateTransformations: IdentityTransformation, Translation
import GeometryTypes
const gt = GeometryTypes
using FileIO
using MeshIO
using Base.Test
import StaticArrays: SVector

const mesh_dir = joinpath(dirname(@__FILE__), "meshes")

@testset "mesh distance" begin
    mesh = load(joinpath(mesh_dir, "base_link.obj"))
    for x in linspace(0.1, 1, 10)
        for y in linspace(-1, 1, 10)
            for z in linspace(-1, 1, 10)
                point = SVector(x, y, z)
                gjk_dist = gjk(mesh, point)
                simple_dist = ReferenceDistance.signed_distance(mesh, point)
                @test isapprox(gjk_dist.signed_distance, simple_dist, atol=2e-4)
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
    for x in linspace(-width, width, 21)
        for y in linspace(-width, width, 21)
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

# @testset "benchmarks" begin
#     include("../perf/runbenchmarks.jl")
# end
