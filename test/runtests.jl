using EnhancedGJK
import EnhancedGJK: projection_weights, projection_weights_reference
using Base.Test
import StaticArrays: SVector

@testset "2d weights" begin
    simplex = SVector{3}(SVector{2, Float64}[[1., -1], [1., 1], [-1.5, 0]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.3, 0.3, 0.4])
    @test isapprox(weights, projection_weights_reference(simplex))

    simplex = SVector{3}(SVector{2, Float64}[[1., -1], [1., 1], [1.5, 0]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.5, 0.5, 0.0])
    @test isapprox(weights, projection_weights_reference(simplex))

    simplex = SVector{3}(SVector{2, Float64}[[1., -1], [1., 1], [0.5, 0]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.0, 0.0, 1.0])
    @test isapprox(weights, projection_weights_reference(simplex))
end

@testset "3d weights" begin
    simplex = SVector{4}(SVector{3, Float64}[[1., -1, 0], [1., 1, 0], [0.4, 0, -0.1], [0.5, 0, 1]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.0, 0.0,0.942623, 0.057377], atol=1e-6)
    @test isapprox(weights, projection_weights_reference(simplex))

    simplex = SVector{4}(SVector{3, Float64}[[1., -1, 0], [1., 1, 0], [1.4, 0, -0.1], [1.5, 0, 1]])
    weights = projection_weights(simplex)
    @test isapprox(weights, [0.5, 0.5, 0.0, 0.0])
    @test isapprox(weights, projection_weights_reference(simplex))
end

@testset "random simplex" begin
    srand(1)
    for i in 1:100
        simplex = SVector{4}([rand(SVector{3, Float64}) for i in 1:4])
        weights = projection_weights(simplex)
        @test all(weights .>= 0)
        @test isapprox(sum(weights), 1)
        @test isapprox(weights, projection_weights_reference(simplex))
    end
end

@testset "benchmarks" begin
    include("../perf/runbenchmarks.jl")
end
