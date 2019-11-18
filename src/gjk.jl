struct Difference{PA, PB}
    a::PA
    b::PB
end

Base.:*(n::Number, d::Difference) = Difference(*(n, d.a), *(n, d.b))
Base.:+(d1::Difference, d2::Difference) = Difference(d1.a + d2.a, d1.b + d2.b)
any_inside(d::Difference) = Difference(any_inside(d.a), any_inside(d.b))

struct CollisionCache{GeomA, GeomB, M, D <: Difference}
    bodyA::GeomA
    bodyB::GeomB
    simplex_points::MVector{M, D}
end

function CollisionCache(geomA, geomB)
    N = dimension(typeof(geomA))
    @assert dimension(typeof(geomA)) == dimension(typeof(geomB))
    CollisionCache(N, geomA, geomB)
end

function CollisionCache(::Val{N}, geomA, geomB) where N
    interior_point = any_inside(Difference(geomA, geomB))
    simplex = MVector(ntuple(i -> interior_point, Val(N + 1)))
    CollisionCache(geomA, geomB, simplex)
end

function reset!(cache::CollisionCache)
    interior_point = any_inside(Difference(cache.bodyA, cache.bodyB))
    simplex = cache.simplex_points
    @inbounds for i in eachindex(simplex)
        simplex[i] = interior_point
    end
    cache
end

dimension(::Type{CollisionCache{G1, G2, M, D}}) where {G1,G2,M,D} = dimension(G1)

function support_vector_max(geometry::gt.GeometryPrimitive, direction, initial_guess::Tagged)
    best_pt, score = gt.support_vector_max(geometry, gt.Vec(direction))
    Tagged(SVector(best_pt))
end

function support_vector_max(pt::gt.Vec{N, T}, direction, initial_guess::Tagged) where {N,T}
    Tagged(SVector(pt))
end

function support_vector_max(pt::gt.Point{N, T}, direction, initial_guess::Tagged) where {N,T}
    Tagged(SVector(pt))
end

function support_vector_max(simplex::Union{gt.AbstractSimplex, gt.AbstractFlexibleGeometry}, direction, initial_guess::Tagged)
    best_pt, score = gt.support_vector_max(simplex, gt.Vec(direction))
    Tagged(SVector(best_pt))
end

function support_vector_max(mesh::gt.HomogenousMesh{gt.Point{N, T}}, direction, initial_guess::Tagged) where {N,T}
    best_arg, best_value = gt.argmax(x-> dot(SVector(x), direction), gt.vertices(mesh))
    best_vec = SVector(best_arg)
    Tagged(best_vec)
end

any_inside(pt::SVector) = Tagged(pt)
support_vector_max(pt::SVector, direction, initial_guess::Tagged) = Tagged(pt)

function transform_simplex(cache::CollisionCache, poseA, poseB)
    transform_simplex(dimension(typeof(cache)), cache, poseA, poseB)
end

@generated function transform_simplex(::Val{N}, cache::CollisionCache, poseA, poseB) where {N}
    transform_simplex_impl(N, cache, poseA, poseB)
end

function transform_simplex_impl(N, cache, poseA, poseB)
    Expr(:call, :(SVector),
        [:((poseA(value(cache.simplex_points[$i].a)) -
            poseB(value(cache.simplex_points[$i].b)))) for i in 1:(N + 1)]...)
end

# Note: it looks like this can be replaced with transpose(weights) * points in Julia 1.3 (before that, it's a lot slower)
@generated function linear_combination(weights::StaticVector{N}, points::StaticVector{N}) where {N}
    expr = :(weights[1] * points[1])
    for i = 2 : N
        expr = :($expr + weights[$i] * points[$i])
    end
    return quote
        Base.@_inline_meta
        $expr
    end
end

struct GJKResult{M, N, T}
    simplex::SVector{M, SVector{N, T}}
    weights::SVector{M, T}
    in_collision::Bool
    closest_point_in_body::Difference{SVector{N, T}, SVector{N, T}}
    iterations::Int
end

function Base.getproperty(result::GJKResult, sym::Symbol)
    if sym === :signed_distance
        @warn "The `signed_distance` field was removed. Please use the `separation_distance` and `simplex_penetration_distance` functions."
        if result.in_collision
            return -simplex_penetration_distance(result)
        else
            return separation_distance(result)
        end
    else
        return Core.getfield(result, sym)
    end
end

closest_point_in_world(result::GJKResult) = linear_combination(result.simplex, result.weights)
closest_point_in_body(result::GJKResult) = result.closest_point_in_body # just for symmetry with the world function
separation_distance(result::GJKResult) = (@assert !result.in_collision; norm(closest_point_in_world(result)))
simplex_penetration_distance(result::GJKResult) = penetration_distance(result.simplex)

function gjk!(cache::CollisionCache,
              poseA::Transformation,
              poseB::Transformation,
              max_iter=100,
              atol=1e-6)
    rotAinv = transform_deriv(inv(poseA), 0)
    rotBinv = transform_deriv(inv(poseB), 0)
    simplex = transform_simplex(cache, poseA, poseB)
    iter = 1

    while true
        weights = projection_weights(simplex)
        min_weight, index_to_replace = findmin(weights)
        best_point = linear_combination(weights, simplex)
        if min_weight > 0 || best_point == zero(best_point)
            # `simplex` contains the origin
            # Nominally in collision; but check for numerical issues
            separation_squared = best_point ⋅ best_point
            @assert separation_squared < sqrt(1_000_000 * eps(typeof(separation_squared)))
            closest_point_in_body = linear_combination(weights, cache.simplex_points)
            return GJKResult(simplex, weights, true, closest_point_in_body, iter)
        end

        direction = -best_point
        direction_in_A = rotAinv * direction
        direction_in_B = rotBinv * direction

        starting_vertex_index = 1
        starting_vertex = cache.simplex_points[starting_vertex_index]
        starting_vertex_score =
            dot(value(starting_vertex.a), direction_in_A) +
            dot(value(starting_vertex.b), direction_in_B)
        for j in 2:length(cache.simplex_points)
            candidate = cache.simplex_points[j]
            candidate_score =
                dot(value(candidate.a), direction_in_A) +
                dot(value(candidate.b), direction_in_B)
            if candidate_score > starting_vertex_score
                starting_vertex_score = candidate_score
                starting_vertex = candidate
            end
        end

        improved_vertex = Difference(
            support_vector_max(cache.bodyA, direction_in_A, starting_vertex.a),
            support_vector_max(cache.bodyB, -direction_in_B, starting_vertex.b))
        improved_point = poseA(value(improved_vertex.a)) - poseB(value(improved_vertex.b))
        score = dot(improved_point, direction)
        if score <= dot(best_point, direction) + atol || iter >= max_iter
            closest_point_in_body = linear_combination(weights, cache.simplex_points)
            return GJKResult(simplex, weights, false, closest_point_in_body, iter)
        else
            cache.simplex_points[index_to_replace] = improved_vertex
            simplex = setindex(simplex, improved_point, index_to_replace)
            # simplex[index_to_replace] = improved_point
        end
        iter += 1
    end
end

function penetration_distance(simplex)
    best_dist_squared = nothing
    for i in eachindex(simplex)
        face = simplex_face(simplex, i)
        weights = projection_weights(face)
        closest_point = linear_combination(weights, face)
        dist_squared = closest_point ⋅ closest_point
        if best_dist_squared === nothing || dist_squared < best_dist_squared
            best_dist_squared = dist_squared
        end
    end
    return sqrt(best_dist_squared)
end

function gjk(geomA, geomB,
             poseA::Transformation=IdentityTransformation(),
             poseB::Transformation=IdentityTransformation())
    cache = CollisionCache(geomA, geomB)
    gjk!(cache, poseA, poseB)
end
