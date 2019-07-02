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


"""
    normal(face::StaticVector{N, <:StaticVector{N}})

Return any `N`-vector normal to the face spanned by `N` points. The result is not necessarily normalized,
nor is there any guarantee on winding.
"""
function normal(face::StaticVector{N, <:StaticVector{N}}) where N
    Δ = hcat(ntuple(i -> face[i + 1] - face[1], Val(N - 1))...)
    factorization = svd(Δ, full=Val(true))
    factorization.U[:, end]
end

function normal(face::StaticVector{2, <:StaticVector{2}})
    Δ = face[2] - face[1]
    similar_type(Δ)(-Δ[2], Δ[1])
end

function normal(face::StaticVector{3, <:StaticVector{3}})
    Δ1 = face[2] - face[1]
    Δ2 = face[3] - face[1]
    Δ1 × Δ2
end
