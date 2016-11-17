using BenchmarkTools
using EnhancedGJK
import EnhancedGJK: projection_weights
import StaticArrays: SVector

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

    suite["GeometryTypes comparisons"] = BenchmarkGroup()
    let
        group = suite["GeometryTypes comparisons"]
        c1 = gt.Simplex(gt.Vec(0.,0), gt.Vec(0.,1), gt.Vec(1.,0),gt.Vec(1.,1))
        c2 = gt.Simplex(gt.Vec(4.,0.5))
        group["geometrytypes simplex"] = @benchmarkable gt.gjk($c1, $c2)

        cache = CollisionCache(c1, c2)
        group["enhanced simplex"] = @benchmarkable gjk!($cache, IdentityTransformation(), IdentityTransformation())
    end

    tune!(suite)
    results = run(suite, verbose=true)
    showall(results)
end
