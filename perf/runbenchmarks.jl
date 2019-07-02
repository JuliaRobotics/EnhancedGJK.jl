using BenchmarkTools
using EnhancedGJK
import GeometryTypes
using CoordinateTransformations: IdentityTransformation, Translation, AffineMap
using EnhancedGJK: projection_weights, transform_simplex, penetration_distance
using StaticArrays: SVector, SMatrix
using Random

Random.seed!(1)

const gt = GeometryTypes

function comparison_benchmark(c1, c2)
    group = BenchmarkGroup()
    cache = CollisionCache(c1, c2)
    group["geometrytypes"] = @benchmarkable gt.gjk($c1, $c2)
    group["enhanced"] = @benchmarkable gjk!($cache, IdentityTransformation(), IdentityTransformation())
    group["enhanced no cache"] = @benchmarkable gjk($c1, $c2)
    return group
end

let
    suite = BenchmarkGroup()

    suite["3d_simplex_1"] = @benchmarkable projection_weights(simplex) setup=(
        simplex = SVector{4}(SVector{3, Float64}[[1., -1, 0], [1., 1, 0], [0.4, 0, -0.1], [0.5, 0, 1]])
    )

    suite["3d_simplex_2"] = @benchmarkable projection_weights(simplex) setup=(
        simplex = SVector{4}(SVector{3, Float64}[[1., -1, 0], [1., 1, 0], [1.4, 0, -0.1], [1.5, 0, 1]])
    )

    suite["2d_simplex_1"] = @benchmarkable projection_weights(simplex) setup=(
        simplex = SVector{3}(SVector{2, Float64}[[1., -1], [1., 1], [0.5, 0]])
    )

    suite["2d_simplex_2"] = @benchmarkable projection_weights(simplex) setup=(
        simplex = SVector{3}(SVector{2, Float64}[[1., -1], [1., 1], [1.5, 0]])
    )

    comparisons = suite["GeometryTypes comparisons"] = BenchmarkGroup()
    comparisons["2d_simplex_simplex"] = let
        c1 = gt.Simplex(gt.Vec(0.,0), gt.Vec(0.,1), gt.Vec(1.,0),gt.Vec(1.,1))
        c2 = gt.Simplex(gt.Vec(4.,0.5))
        comparison_benchmark(c1, c2)
    end
    comparisons["3d_point_box"] = let
        c1 = gt.HyperRectangle(rand(gt.Vec{3}), rand(gt.Vec{3}))
        c2 = rand(gt.Vec{3})
        comparison_benchmark(c1, c2)
    end

    suite["transform_simplex"] = @benchmarkable transform_simplex(cache, poseA, poseB) setup = begin
        c1 = gt.HyperRectangle(rand(gt.Vec{3}), rand(gt.Vec{3}))
        c2 = gt.Simplex(rand(gt.Vec{3}))
        cache = CollisionCache(c1, c2)
        poseA = AffineMap(rand(SMatrix{3, 3}), rand(SVector{3}))
        poseB = AffineMap(rand(SMatrix{3, 3}), rand(SVector{3}))
        gjk!(cache, poseA, poseB)
    end

    suite["penetration_distance"] = @benchmarkable penetration_distance(simplex) setup = begin
        c1 = gt.HyperRectangle(rand(gt.Vec{3}), rand(gt.Vec{3}))
        c2 = rand(gt.Vec{3})
        result = gjk(c1, c2)
        simplex = result.simplex
    end

    tune!(suite)
    results = run(suite)
    show(IOContext(stdout, :compact => false), results)
end
