using EnhancedGJK
import EnhancedGJK: projection_weights
using Base.Test
import StaticArrays: SVector

@testset "2d weights" begin
    simplex = SVector{3}(SVector{2, Float64}[[1., -1], [1., 1], [-1.5, 0]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.3, 0.3, 0.4])

    simplex = SVector{3}(SVector{2, Float64}[[1., -1], [1., 1], [1.5, 0]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.5, 0.5, 0.0])

    simplex = SVector{3}(SVector{2, Float64}[[1., -1], [1., 1], [0.5, 0]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.0, 0.0, 1.0])
end

@testset "3d weights" begin
    simplex = SVector{4}(SVector{3, Float64}[[1., -1, 0], [1., 1, 0], [0.4, 0, -0.1], [0.5, 0, 1]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.0, 0.0,0.942623, 0.057377], atol=1e-6)

    simplex = SVector{4}(SVector{3, Float64}[[1., -1, 0], [1., 1, 0], [1.4, 0, -0.1], [1.5, 0, 1]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.5, 0.5, 0.0, 0.0])
end

@testset "benchmarks" begin
    include("../perf/runbenchmarks.jl")
end
