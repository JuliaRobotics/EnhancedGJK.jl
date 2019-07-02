num_johnson_subsets(simplex_length::Integer) = 2^simplex_length - 1

"""
Compute all subsets of the points in the simplex in a reliable order. The
order is arbitrary, but was chosen to match the implemenation in
S. Cameron, “Enhancing GJK: computing minimum and penetration distances between
convex polyhedra,”. Specifically, subset i contains point j iff the binary
representation of i has a one at the jth bit.
"""
function johnson_subsets(simplex_length::Integer)
    num_subsets = num_johnson_subsets(simplex_length)
    subsets = falses(simplex_length, num_subsets)
    for i in 1:num_subsets
        for j in 1:simplex_length
            subsets[j, i] = (i >> (j - 1)) & 1
        end
    end
    subsets
end

"""
    weights = projection_weights(simplex)

This function implements Johnson's distance subalgorithm, as described in
E. G. Gilbert, D. W. Johnson, and S. S. Keerthi, “A fast procedure for
computing the distance between complex objects in three-dimensional space,”
1988. Given a simplex (a length N+1 vector of points of dimension N), it
returns weights such that dot(weights, simplex) yields the point in the convex
hull of the simplex which is closest to the origin.

This is the critical loop of the GJK algorithm, so it has been heavily optimized
to precompute, inline, and unroll as much of the algorithm as possible. For a
more readable (and much slower) implementation, see
projection_weights_reference()
"""
@generated function projection_weights(simplex::SVector{M, SVector{N, T}}; atol=sqrt(eps(T))) where {M,N,T}
    simplex_length = M
    num_subsets = num_johnson_subsets(M)
    subsets = johnson_subsets(M)
    complements = .!subsets

    expr = quote
        deltas = SVector(tuple($([zero(SVector{simplex_length, T}) for i = 1 : num_subsets]...)))
    end

    # Set the weight of every singleton subset to 1.
    for i in 1:simplex_length
        elements = [:(zero($T)) for j in 1:simplex_length]
        elements[i] = :(one($T))
        arg = :(SVector{$simplex_length, $T}(tuple($(elements...))))

        push!(expr.args, quote
            deltas = setindex(deltas, $(arg), $(2^(i - 1)))
        end)
    end

    for s in 1:(num_subsets - 1)
        k = first((1:M)[subsets[:,s]])
        push!(expr.args, quote
            viable = true
        end)
        for j in (1:M)[complements[:,s]]
            s2 = s + (1 << (j - 1))
            push!(expr.args, quote
                d = deltas[$s2][$j]
            end)
            for i in (1:M)[subsets[:,s]]
                push!(expr.args, quote
                    d += deltas[$s][$i] *
                    (dot(simplex[$i], simplex[$k]) - dot(simplex[$i], simplex[$j]))
                end)
            end
            push!(expr.args, quote
                deltas = setindex(deltas, setindex(deltas[$s2], d, $j), $s2)
            end)
            push!(expr.args, quote
                if d > atol # See #17
                    viable = false
                end
            end)
        end
        push!(expr.args, quote
            viable = viable && all($(Expr(:tuple,
            [:(deltas[$s][$i] >= 0) for i in (1:M)[subsets[:,s]]]...)))
            if viable
                return deltas[$s] ./ sum(deltas[$s])
            end
        end)
    end
    push!(expr.args, quote
        return deltas[$num_subsets] ./ sum(deltas[$num_subsets])
    end)
    return expr
end

function projection_weights_reference(simplex::SVector{M, SVector{N, T}}) where {M,N,T}
    num_subsets = 2^M - 1
    complements = falses(M, num_subsets)
    subsets = falses(M, num_subsets)
    for i in 1:num_subsets
        for j in 1:M
            subsets[j, i] = (i >> (j - 1)) & 1
            complements[j, i] = !subsets[j, i]
        end
    end

    deltas = zeros(MMatrix{M, num_subsets, T})
    viable = ones(MVector{num_subsets, Bool})

    for i in 1:M
        deltas[i, 2^(i - 1)] = 1
    end

    for s in 1:(num_subsets - 1)
        k = first((1:M)[subsets[:,s]])
        for j in (1:M)[complements[:,s]]
            s2 = s + (1 << (j - 1))
            delta_j_s2 = 0
            for i in (1:M)[subsets[:,s]]
                delta_j_s2 += deltas[i, s] * (dot(simplex[i], simplex[k]) - dot(simplex[i], simplex[j]))
            end
            if delta_j_s2 > 0
                viable[s] = false
            elseif delta_j_s2 < 0
                viable[s2] = false
            end
            deltas[j, s2] = delta_j_s2
        end
        if viable[s]
            return deltas[:,s] ./ sum(deltas[:,s])
        end
    end
    return deltas[:,end] ./ sum(deltas[:,end])
end
