typealias Simplex{M, N, T} Union{MVector{M, SVector{N, T}}, SVector{M, SVector{N, T}}}
typealias SimplexType{M, N, T} Union{Type{MVector{M, SVector{N, T}}}, Type{SVector{M, SVector{N, T}}}}

@pure dimension{M, N, T}(::SimplexType{M, N, T}) = Val{N}
@pure scalartype{M, N, T}(::SimplexType{M, N, T}) = T

any_inside{M, N, T}(simplex::Simplex{M, N, T}) = Tagged(simplex[1])

function support_vector_max{M, N, T}(simplex::Simplex{M, N, T}, direction, initial_guess::Tagged)
    best_pt, score = gt.argmax(p -> dot(p, direction), simplex)
    Tagged(best_pt)
end
