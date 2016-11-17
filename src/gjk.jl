immutable Difference{PA, PB}
    a::PA
    b::PB
end

*(n::Number, d::Difference) = Difference(*(n, d.a), *(n, d.b))
+(d1::Difference, d2::Difference) = Difference(d1.a + d2.a, d1.b + d2.b)

immutable CollisionCache{GeomA, GeomB, M, D <: Difference}
    bodyA::GeomA
    bodyB::GeomB
    simplex_points::MVector{M, D}
end

function CollisionCache(geomA, geomB)
    N = dimension(typeof(geomA))
    @assert dimension(typeof(geomA)) == dimension(typeof(geomB))
    CollisionCache(N, geomA, geomB)
end

@generated function initial_simplex{N}(dim::Type{Val{N}}, geomA, geomB)
    return quote
        interior_point = Difference(any_inside(geomA), any_inside(geomB))
        simplex = $(Expr(:call, :(MVector), [:(interior_point) for i in 1:(N + 1)]...))
    end
end

function CollisionCache{N}(::Type{Val{N}}, geomA, geomB)
    simplex = initial_simplex(Val{N}, geomA, geomB)
    CollisionCache(geomA, geomB, simplex)
end

dimension{G1, G2, M, D}(::Type{CollisionCache{G1, G2, M, D}}) = dimension(G1)

function support_vector_max(geometry::gt.GeometryPrimitive, direction, initial_guess::Tagged)
    best_pt, score = gt.support_vector_max(geometry, gtvec(direction))
    Tagged(svector(best_pt))
end

function support_vector_max{N, T}(pt::gt.Vec{N, T}, direction, initial_guess::Tagged)
    Tagged(svector(pt))
end

function support_vector_max(simplex::Union{gt.AbstractSimplex, gt.AbstractFlexibleGeometry}, direction, initial_guess::Tagged)
    best_pt, score = gt.support_vector_max(simplex, gtvec(direction))
    Tagged(svector(best_pt))
end

function support_vector_max{N, T}(mesh::gt.HomogenousMesh{gt.Point{N, T}}, direction, initial_guess::Tagged)
    best_arg, best_value = gt.argmax(x-> dot(svector(x), direction), gt.vertices(mesh))
    best_vec = svector(best_arg)
    Tagged(best_vec)
end

any_inside(pt::SVector) = Tagged(pt)
support_vector_max(pt::SVector, direction, initial_guess::Tagged) = Tagged(pt)

function transform_simplex(cache::CollisionCache, poseA, poseB)
    transform_simplex(dimension(typeof(cache)), cache, poseA, poseB)
end

@generated function transform_simplex{N}(::Type{Val{N}}, cache::CollisionCache, poseA, poseB)
    transform_simplex_impl(N, cache, poseA, poseB)
end

function transform_simplex_impl(N, cache, poseA, poseB)
    Expr(:call, :(SVector),
        [:((poseA(value(cache.simplex_points[$i].a)) -
            poseB(value(cache.simplex_points[$i].b)))) for i in 1:(N + 1)]...)
end

immutable GJKResult{M, N, T}
    simplex::SVector{M, SVector{N, T}}
    closest_point_in_body::Difference{SVector{N, T}, SVector{N, T}}
    signed_distance::T
end

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
        if min_weight > 0
            # in collision
            return GJKResult(
                simplex,
                dot(weights, cache.simplex_points),
                penetration_distance(simplex)
            )
        end
        best_point = dot(weights, simplex)

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
            return GJKResult(
                simplex,
                dot(weights, cache.simplex_points),
                norm(best_point)
            )
        else
            cache.simplex_points[index_to_replace] = improved_vertex
            simplex = setindex(simplex, improved_point, index_to_replace)
            # simplex[index_to_replace] = improved_point
        end
        iter += 1
    end
end

function penetration_distance(simplex)
    _, penetration_distance = gt.argmax(1:length(simplex)) do i
        face = simplex_face(simplex, i)
        weights = projection_weights(face)
        closest_point = dot(weights, face)
        distance_to_face = norm(closest_point)
        -distance_to_face
    end
    return penetration_distance
end

function gjk(geomA, geomB,
             poseA::Transformation=IdentityTransformation(),
             poseB::Transformation=IdentityTransformation())
    cache = CollisionCache(geomA, geomB)
    gjk!(cache, poseA, poseB)
end
