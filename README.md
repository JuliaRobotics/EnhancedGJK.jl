# EnhancedGJK

[![Build Status](https://travis-ci.org/rdeits/EnhancedGJK.jl.svg?branch=master)](https://travis-ci.org/rdeits/EnhancedGJK.jl)
[![codecov.io](http://codecov.io/github/rdeits/EnhancedGJK.jl/coverage.svg?branch=master)](http://codecov.io/github/rdeits/EnhancedGJK.jl?branch=master)

This package contains a pure-Julia implementation of the enhanced version of the Gilbert, Johnson, and Keerthi algorithm for computing distance between convex bodies. The algorithm is described in detail by Stephen Cameron in [1].

## Why Julia?

GJK implementations are numerous and well-tested, but a pure-Julia implementation may have benefits that other languages cannot offer. This implementation of GJK is entirely agnostic to the data types which describe both the geometries and their positions in space. This means that, for example, gradients of distances can easily be computed using the automatic differentiation provided by ForwardDiff.jl's DualNumber type. But there may be other applications as well, such as geometries or transformations with rational or variable-precision arithmetic. A pure-Julia implementation makes it easy to experiment with new data types without sacrificing performance.

# Usage

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
# Reusing the same cache will make this computation faster for complex
# geometries if the change in transformation is small.
result = gjk!(cache, Translation(SVector(0.1, 0, 0)), IdentityTransformation())

@show result.signed_distance
```


# References

[1] S. Cameron, “Enhancing GJK: computing minimum and penetration distances between convex polyhedra,” in , 1997 IEEE International Conference on Robotics and Automation, 1997. Proceedings, 1997, vol. 4, pp. 3112–3117 vol.4.
