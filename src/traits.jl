dimension(::Type{gt.AbstractSimplex{S, gt.Vec{N, T}}}) where {N,T,S} = Val{N}

dimension(::Type{M}) where {M <: gt.AbstractMesh} = dimension(gt.vertextype(M))

dimension(::Type{gt.Vec{N, T}}) where {N,T} = Val{N}

dimension(::Type{gt.Point{N, T}}) where {N,T} = Val{N}

dimension(::Type{gt.Simplex{M, T}}) where {M,T} = dimension(T)

dimension(::Type{SVector{N, T}}) where {N,T} = Val{N}

dimension(::Type{gt.FlexibleConvexHull{T}}) where {T} = dimension(T)

dimension(::Type{gt.HyperRectangle{N, T}}) where {N,T} = Val{N}
