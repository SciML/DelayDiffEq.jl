"""
    Discontinuity(t, order::Int)

Object of discontinuity of order `order` at time `t`, i.e. discontinuity of
`order`th derivative at time `t`.
"""
struct Discontinuity{tType}
    t::tType
    order::Int
end

# ordering of discontinuities
Base.:<(a::Discontinuity, b::Discontinuity) = a.t < b.t || (a.t == b.t && a.order < b.order)
Base.isless(a::Discontinuity, b::Discontinuity) = isless(a.t, b.t) || (isequal(a.t, b.t) && a.order < b.order)
Base.isequal(a::Discontinuity, b::Discontinuity) = isequal(a.t, b.t) && a.order == b.order
Base.:(==)(a::Discontinuity, b::Discontinuity) = a.t == b.t && a.order == b.order

# ordering with numbers
Base.:<(a::Discontinuity, b::Number) = a.t < b
Base.:<(a::Number, b::Discontinuity) = a < b.t
Base.isless(a::Discontinuity, b::Number) = isless(a.t, b)
Base.isless(a::Number, b::Discontinuity) = isless(a, b.t)
Base.:(==)(a::Discontinuity, b::Number) = a.t == b
Base.:(==)(a::Number, b::Discontinuity) = a == b.t
Base.isequal(a::Discontinuity, b::Number) = isequal(a.t, b)
Base.isequal(a::Number, b::Discontinuity) = isequal(a, b.t)

# multiplication with numbers
Base.:*(a::Discontinuity, b::Number) = a.t * b
Base.:*(a::Number, b::Discontinuity) = a * b.t

# conversion to numbers
Base.convert(::Type{T}, d::Discontinuity) where {T<:Number} = convert(T, d.t)
Base.convert(::Type{T}, d::Discontinuity{T}) where {T<:Number} = d.t

# display
Base.show(io::IO, d::Discontinuity) = print(io, d.t, " (order ", d.order, ")")
Base.show(io::IO, ::MIME"text/plain", d::Discontinuity{T}) where {T} =
    print(io, "Discontinuity{", T, "}:\n   ", d)
