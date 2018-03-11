@pure dimension(::Type{gt.AbstractSimplex{S, gt.Vec{N, T}}}) where {N,T,S} = Val{N}

@pure dimension(::Type{M}) where {M <: gt.AbstractMesh} = dimension(gt.vertextype(M))

@pure dimension(::Type{gt.Vec{N, T}}) where {N,T} = Val{N}

@pure dimension(::Type{gt.Point{N, T}}) where {N,T} = Val{N}

@pure dimension(::Type{gt.Simplex{M, T}}) where {M,T} = dimension(T)

@pure dimension(::Type{SVector{N, T}}) where {N,T} = Val{N}

@pure dimension(::Type{gt.FlexibleConvexHull{T}}) where {T} = dimension(T)

@pure dimension(::Type{gt.HyperRectangle{N, T}}) where {N,T} = Val{N}
