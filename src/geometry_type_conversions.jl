# A few hack-y conversions to and from the FixedSizeArrays used
# within GeometryTypes. These can probably all go away when
# GeometryTypes switches to using StaticArrays internally.

@generated function svector(v::Union{gt.Vec{N, T}, gt.Point{N, T}}) where {N,T}
    svector_impl(v)
end

function svector_impl(::Union{Type{gt.Vec{N, T}}, Type{gt.Point{N, T}}}) where {N,T}
    Expr(:call, :(SVector{$N, $T}), [:(v[$i]) for i in 1:N]...)
end

@generated function gtvec(v::SVector{N, T}) where {N,T}
    Expr(:call, :(gt.Vec), [:(v[$i]) for i in 1:N]...)
end
