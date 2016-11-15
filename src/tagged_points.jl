immutable Tagged{P, T}
    point::P
    tag::T
end

Tagged{P}(point::P) = Tagged(point, nothing)

*(n::Number, t::Tagged) = n * value(t)

value(t::Tagged) = t.point

function any_inside(geometry)
    point = gt.any_inside(geometry)
    Tagged(svector(point))
end

function any_inside{N, T}(mesh::gt.AbstractMesh{gt.Point{N, T}})
    Tagged(svector(first(gt.vertices(mesh))))
end

# function any_inside(m::gt.MinkowskiDifference)
#     t1 = any_inside(m.c1)
#     t2 = any_inside(m.c2)
#     Difference(t1, t2)
# end
