@generated function projection_weights{M, N, T}(simplex::SVector{M, SVector{N, T}})
    projection_weights_impl(simplex)
end

num_johnson_subsets(simplex_length::Integer) = 2^simplex_length - 1

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

function initial_deltas_impl(simplex_length, T)
    num_subsets = num_johnson_subsets(simplex_length)

    expr = quote
        deltas = zeros(SVector{$num_subsets, SVector{$simplex_length, $T}})
    end

    # Set the weight of every singleton subset to 1.
    for i in 1:simplex_length
        elements = [:(zero($T)) for j in 1:simplex_length]
        elements[i] = :(one($T))
        arg = Expr(:call, :(SVector{$simplex_length, $T}), elements...)

        push!(expr.args, quote
            deltas = setindex(deltas, $(arg), $(2^(i - 1)))
        end)
    end
    expr
end

"""
This is the function which actually computes the weight values.
It returns an *expression* to compute that result, with all loops
unrolled and all indices pre-computed based on the size and dimension
of the simplex.
"""
function projection_weights_impl{M, N, T}(::Type{SVector{M, SVector{N, T}}})
    num_subsets = num_johnson_subsets(M)
    subsets = johnson_subsets(M)
    complements = !subsets

    expr = quote
        deltas = $(initial_deltas_impl(M, T))
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
                if d > 0
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

function projection_weights_reference{M, N, T}(simplex::SVector{M, SVector{N, T}})
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
