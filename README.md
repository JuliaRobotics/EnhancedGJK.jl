# EnhancedGJK

[![Build Status](https://travis-ci.org/JuliaRobotics/EnhancedGJK.jl.svg?branch=master)](https://travis-ci.org/JuliaRobotics/EnhancedGJK.jl)
[![codecov.io](http://codecov.io/github/JuliaRobotics/EnhancedGJK.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaRobotics/EnhancedGJK.jl?branch=master)

This package contains a pure-Julia implementation of the enhanced version of the Gilbert, Johnson, and Keerthi algorithm for computing distance between convex bodies. The algorithm is described in detail by Stephen Cameron in [1].

## Why Julia?

GJK implementations are numerous and well-tested, but a pure-Julia implementation may have benefits that other languages cannot offer. This implementation of GJK is entirely agnostic to the data types which describe both the geometries and their positions in space. This means that, for example, gradients of distances can easily be computed using the automatic differentiation provided by ForwardDiff.jl's DualNumber type. But there may be other applications as well, such as geometries or transformations with rational or variable-precision arithmetic. A pure-Julia implementation makes it easy to experiment with new data types without sacrificing performance.

# Usage

The easiest way to use this package is the `gjk()` function. `gjk()` takes two geometries and, optionally, Transformation types specifying the pose of each geometry:

```julia
using EnhancedGJK
import GeometryTypes: HyperRectangle, HyperSphere, Vec, Point

c1 = HyperRectangle(Vec(0., 0, 0), Vec(1., 1, 1))
c2 = HyperRectangle(Vec(3., 0, 0), Vec(1., 1, 1))
result = gjk(c1, c2)
```

The return type of `gjk()` is a `GJKResult`, from which you can extract the signed distance between the two bodies:

```julia
julia> @show result.signed_distance
result.signed_distance = 2.0
```

You can also access the closest point in each body to the other:

```julia
julia> result.closest_point_in_body.a
3-element StaticArrays.SVector{3,Float64}:
 1.0
 0.0
 0.0

julia> result.closest_point_in_body.b
3-element StaticArrays.SVector{3,Float64}:
 3.0
 0.0
 0.0
```

## Going Faster

When simulating physics, we often want to compute the distance between two bodies over and over while those bodies move slightly. In that case, we can cache some of the intermediate results to make each distance computation faster and free of memory allocations:

```julia
using EnhancedGJK
import StaticArrays: SVector
import CoordinateTransformations: IdentityTransformation, Translation

# Construct two geometries: a simplex and a single point:
simplex = SVector{3}(SVector{2, Float64}[[1., 0], [2., 0], [1., 1]])
pt = SVector(0., 0)

# The CollisionCache stores both geometries and also remembers
# information about the GJK simplex used to check for collisions
# between them. Using the same cache later will make subsequent
# computations faster.
cache = CollisionCache(simplex, pt);

# Run the GJK algorithm. Each geometry can also be given a
# transformation to describe its position and orientation in the
# world frame.
result = gjk!(cache, IdentityTransformation(), IdentityTransformation())

# result.signed_distance will be > 0 if the objects are not in contact
# and <= 0 if they are in collision.
@show result.signed_distance

# We can perturb one of the geometries by changing its transformation.
# Reusing the same cache will make this computation faster, expecially
# for complex geometries when the change in transformation is small.
result = gjk!(cache, Translation(SVector(0.1, 0)), IdentityTransformation())

@show result.signed_distance
```

## Meshes

`gjk()` and `gjk!()` support meshes, represented as GeometryTypes.jl HomogenousMesh objects:

```julia
using MeshIO
using FileIO
mesh = load("test/meshes/r_foot_chull.obj")
result = gjk(mesh, mesh, IdentityTransformation(), Translation(SVector(5., 0, 0)))
@show result.signed_distance
```

Note that this package *does not* check if the mesh is convex. Non-convex meshes may produce incorrect distance measurements.

GJK can be run even faster on complex meshes by pre-computing the neighbors of each vertex in the mesh. The `NeighborMesh` type handles this for you:

```julia
neighbormesh = NeighborMesh(mesh)
result = gjk(neighbormesh, neighbormesh, IdentityTransformation(), Translation(SVector(5., 0, 0)))
```

This pre-computation of mesh vertex neighbors is the "enhanced" part of Enhanced GJK.


# References

[1] S. Cameron, “Enhancing GJK: computing minimum and penetration distances between convex polyhedra,” in , 1997 IEEE International Conference on Robotics and Automation, 1997. Proceedings, 1997, vol. 4, pp. 3112–3117 vol.4.
