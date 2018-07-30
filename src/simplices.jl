const Simplex{M, N, T} = Union{MVector{M, SVector{N, T}}, SVector{M, SVector{N, T}}}
const SimplexType{M, N, T} = Union{Type{MVector{M, SVector{N, T}}}, Type{SVector{M, SVector{N, T}}}}

dimension(::SimplexType{M, N, T}) where {M,N,T} = Val(N)

any_inside(simplex::Simplex{M, N, T}) where {M,N,T} = Tagged(simplex[1])

function support_vector_max(simplex::Simplex{M, N, T}, direction, initial_guess::Tagged) where {M,N,T}
    best_pt, score = gt.argmax(p -> dot(p, direction), simplex)
    Tagged(best_pt)
end

@generated function simplex_face(simplex::Simplex{M, N, T}, i::Integer) where {N,M,T}
    simplex_face_impl(simplex, i)
end

function simplex_face_impl(simplex::SimplexType{M, N, T}, i) where {N,M,T}
    Expr(:call, :(SVector),
        Expr(:tuple, [:(i > $j ? simplex[$j] : simplex[$(j+1)]) for j in 1:(M - 1)]...))
end
