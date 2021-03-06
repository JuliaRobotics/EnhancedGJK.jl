using Test
using Statistics
using LinearAlgebra
using Random

using EnhancedGJK
using EnhancedGJK: projection_weights, projection_weights_reference, reset!, normal
using CoordinateTransformations: IdentityTransformation, Translation
using StaticArrays: SVector
import GeometryTypes
const gt = GeometryTypes
using FileIO
using ForwardDiff

const mesh_dir = joinpath(dirname(@__FILE__), "meshes")

@testset "Issue #17" begin
    # Note: when this test was written, it mattered whether it was the first
    # test or not.
    geomA = gt.FlexibleConvexHull(gt.Point{2,Float64}[
        gt.Point(1.125, 2.0),
        gt.Point(0.125, 2.0),
        gt.Point(0.125, 2.8),
        gt.Point(0.325, 3.0),
        gt.Point(0.925, 3.0),
        gt.Point(1.125, 2.8)]
    )
    geomB = gt.FlexibleConvexHull(gt.Point{2,Float64}[
        gt.Point(-0.025, 1.2321067811865476),
        gt.Point(-0.025, 1.2821067811865474),
        gt.Point(0.3267766952966369, 1.2821067811865474),
        gt.Point(0.3267766952966369, 1.2321067811865476)]
    )

    cache = CollisionCache(geomA, geomB)
    poseA = poseB = IdentityTransformation()
    simplex = EnhancedGJK.transform_simplex(cache, poseA, poseB)
    @test isapprox(projection_weights(simplex), projection_weights_reference(simplex))

    result = gjk(geomA, geomB)
    @test !result.in_collision
end

@testset "reference distance" begin
    mesh = load(joinpath(mesh_dir, "base_link.obj"))
    for x in range(0.05, stop=1, length=10)
        for y in range(-1, stop=1, length=10)
            for z in range(-1, stop=1, length=10)
                point = SVector(x, y, z)
                result = gjk(mesh, point)
                simple_dist = ReferenceDistance.signed_distance(mesh, point)
                if result.in_collision
                    @test isapprox(simplex_penetration_distance(result), -simple_dist, atol=1e-3)
                else
                    @test isapprox(separation_distance(result), simple_dist, atol=1e-3)
                end
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
            let result = gjk(geometry, point, IdentityTransformation(), Translation(SVector(x, y, z)))
                @test separation_distance(result) ≈ 0.05
            end
            z = 0.06
            let result = gjk(geometry, point, IdentityTransformation(), Translation(SVector(x, y, z)))
                @test separation_distance(result) ≈ 0.01 atol=1e-12
            end
            z = 0.05
            let result = gjk(geometry, point, IdentityTransformation(), Translation(SVector(x, y, z)))
                if result.in_collision
                    @test simplex_penetration_distance(result) ≈ 0.0 atol=1e-12
                else
                    @test separation_distance(result) ≈ 0.0 atol=1e-12
                end
            end
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
    @test isapprox(separation_distance(result), 1.0)
    @test isapprox(result.closest_point_in_body.a, [1.0, 0.0])
    @test isapprox(result.closest_point_in_body.b, [0.0, 0.0])
end

@testset "mesh to mesh" begin
    mesh = load(joinpath(mesh_dir, "r_foot_chull.obj"))
    dx = 1.0
    foot_length = 0.172786 + 0.090933

    cache = CollisionCache(mesh, mesh)
    result = gjk!(cache, IdentityTransformation(), Translation(SVector(dx, 0, 0)))
    @test isapprox(separation_distance(result), dx - foot_length, atol=1e-3)

    cache = CollisionCache(mesh, mesh)
    result = gjk!(cache, Translation(SVector(dx, 0, 0)), IdentityTransformation())
    @test isapprox(separation_distance(result), dx - foot_length, atol=1e-3)

    cache = CollisionCache(mesh, mesh)
    expected_penetration = 0.01
    result = gjk!(cache, IdentityTransformation(), Translation(foot_length - expected_penetration, 0, 0))
    # TODO: penetration distance is inconsistent and inaccurate
    @test simplex_penetration_distance(result) > 0
end

@testset "neighbor mesh to mesh" begin
    mesh = NeighborMesh(load(joinpath(mesh_dir, "r_foot_chull.obj")))
    dx = 1.0
    foot_length = 0.172786 + 0.090933

    cache = CollisionCache(mesh, mesh)
    result = gjk!(cache, IdentityTransformation(), Translation(SVector(dx, 0, 0)))
    @test isapprox(separation_distance(result), dx - foot_length, atol=1e-3)

    cache = CollisionCache(mesh, mesh)
    result = gjk!(cache, Translation(SVector(dx, 0, 0)), IdentityTransformation())
    @test isapprox(separation_distance(result), dx - foot_length, atol=1e-3)

    cache = CollisionCache(mesh, mesh)
    expected_penetration = 0.01
    result = gjk!(cache, IdentityTransformation(), Translation(foot_length - expected_penetration, 0, 0))
    # TODO: penetration distance is inconsistent and inaccurate
    @test simplex_penetration_distance(result) > 0
end

@testset "geometry types" begin
    # Adapted from
    # https://github.com/JuliaGeometry/GeometryTypes.jl/blob/master/test/gjk.jl
    @testset "gjk examples" begin
        c1 = gt.Simplex(gt.Vec(-1.))
        c2 = gt.Simplex(gt.Vec(4.))
        cache = CollisionCache(c1, c2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(separation_distance(result), 5.0)

        c1 = gt.Simplex(gt.Vec(-1.,0,0))
        c2 = gt.Simplex(gt.Vec(4.,0,0))
        cache = CollisionCache(c1, c2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(separation_distance(result), 5.0)

        c1 = gt.FlexibleConvexHull([gt.Vec(0.,0), gt.Vec(0.,1), gt.Vec(1.,0),gt.Vec(1.,1)])
        c2 = gt.Simplex(gt.Vec(4.,0.5))
        cache = CollisionCache(c1, c2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(separation_distance(result), 3.0)

        pt1 = gt.Vec(1,2,3.)
        pt2 = gt.Vec(3,4,5.)
        cache = CollisionCache(pt1, pt2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(separation_distance(result), norm(pt1 - pt2))

        pt1 = gt.Point(1,2,3.)
        pt2 = gt.Point(3,4,5.)
        cache = CollisionCache(pt1, pt2)
        result = gjk!(cache, IdentityTransformation(), IdentityTransformation())
        @test isapprox(separation_distance(result), norm(pt1 - pt2))
    end

    @testset "gjk intersecting lines" begin
        c1 = gt.Simplex(gt.Vec(1,1.), gt.Vec(1, 2.))
        @test simplex_penetration_distance(gjk(c1, c1)) == 0.

        c2 = gt.Simplex(gt.Vec(1,1.), gt.Vec(10, 2.))
        @test simplex_penetration_distance(gjk(c1, c2)) == 0.
    end
end

@testset "reset!" begin
    c1 = gt.FlexibleConvexHull([gt.Vec(0.,0), gt.Vec(0.,1), gt.Vec(1.,0),gt.Vec(1.,1)])
    c2 = gt.Simplex(gt.Vec(4.,0.5))
    cache = CollisionCache(c1, c2)
    initial_simplex = deepcopy(cache.simplex_points)
    gjk!(cache, IdentityTransformation(), IdentityTransformation())
    @test cache.simplex_points != initial_simplex
    reset!(cache)
    @test cache.simplex_points == initial_simplex
    allocs = @allocated reset!(cache)
    @test allocs == 0
end

@testset "normal, N = $N" for N = 2 : 4
    rng = MersenneTwister(1)
    for j = 1 : 100
        face = SVector{N}(ntuple(i -> rand(rng, SVector{N}), Val(N)))
        n = normal(face)
        dot_products = [point ⋅ n for point in face]
        @test all(x -> isapprox(x, 0; atol=1e-10), dot_products .- mean(dot_products))
    end
end

@testset "starting simplex is origin" begin
    hr = gt.HyperRectangle{3, Float64}([-1.0, -1.0, -1.0], [2.0, 2.0, 2.0])
    result = gjk(hr, hr)
    @test result.in_collision
end

@testset "Issue #36" begin
    function distance_from_segment(z)
        p1 = GeometryTypes.Point(4.0, -0.5)
        p2 = GeometryTypes.Point(6.0, 0.0)
        l = GeometryTypes.LineSegment(p1, p2)
        p = GeometryTypes.Point(5.0, z)
        result = EnhancedGJK.gjk(l, p)
        return result.in_collision ? 0.0 : separation_distance(result)
    end

    z = -1.0
    deriv_autodiff = ForwardDiff.derivative(distance_from_segment, z)

    # Check the autodiff result against finite difference
    δz = sqrt(eps(Float64))
    deriv_numeric = (distance_from_segment(z + δz) - distance_from_segment(z)) / δz
    @test deriv_autodiff ≈ deriv_numeric
end

@testset "benchmarks" begin
    include("../perf/runbenchmarks.jl")
end
