# EnhancedGJK

[![Build Status](https://travis-ci.org/rdeits/EnhancedGJK.jl.svg?branch=master)](https://travis-ci.org/rdeits/EnhancedGJK.jl)
[![codecov.io](http://codecov.io/github/rdeits/EnhancedGJK.jl/coverage.svg?branch=master)](http://codecov.io/github/rdeits/EnhancedGJK.jl?branch=master)

This package contains a pure-Julia implementation of the enhanced version of the Gilbert, Johnson, and Keerthi algorithm for computing distance between convex bodies. The algorithm is described in detail by Stephen Cameron in [1].

## Why Julia?

GJK implementations are numerous and well-tested, but a pure-Julia implementation may have benefits that other languages cannot offer. This implementation of GJK is entirely agnostic to the data types which describe both the geometries and their positions in space. This means that, for example, gradients of distances can easily be computed using the automatic differentiation provided by ForwardDiff.jl's DualNumber type. But there may be other applications as well, such as geometries or transformations with rational or variable-precision arithmetic. A pure-Julia implementation makes it easy to experiment with new data types without sacrificing performance. 

# References

[1] S. Cameron, “Enhancing GJK: computing minimum and penetration distances between convex polyhedra,” in , 1997 IEEE International Conference on Robotics and Automation, 1997. Proceedings, 1997, vol. 4, pp. 3112–3117 vol.4.
