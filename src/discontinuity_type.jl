"""
    Discontinuity(t, order::Int)

Object of discontinuity of order `order` at time `t`.
"""
struct Discontinuity{tType<:Number} <: Number
    t::tType
    order::Int
end

# sorting
Base.isless(a::Discontinuity, b::Discontinuity) = isless(a.t, b.t)
Base.isless(a::Discontinuity, b::Number) = isless(a.t, b)
Base.isless(a::Number, b::Discontinuity) = isless(a, b.t)

# conversion
Base.convert(::Type{T}, d::Discontinuity{T}) where {T<:Number} = d.t
Base.convert(::Type{T}, d::Discontinuity) where {T<:Number} = convert(T, d.t)
Base.convert(::Type{D}, d::Discontinuity) where {D<:Discontinuity} =
    D(d.t, d.order)

# avoid ambigous methods
Base.convert(::Type{Quantity{T,D,U}}, d::Discontinuity{Quantity{T,D,U}}) where {T,D,U} =
    d.t
Base.convert(::Type{Quantity{T}}, d::Discontinuity) where {T} =
    convert(Quantity{T}, d.t)

# promotion
Base.promote_rule(::Type{Discontinuity{S}}, ::Type{T}) where {S<:Number,T<:Number} =
    promote_type(S, T)
Base.promote_rule(::Type{Discontinuity{S}}, ::Type{Discontinuity{T}}) where
    {S<:Number,T<:Number} =
        Discontinuity{promote_type(S, T)}

# display
Base.show(io::IO, d::Discontinuity) = print(io, d.t, " (order: ", d.order, ")")
Base.show{T}(io::IO, ::MIME"text/plain", d::Discontinuity{T}) =
    print(io, "Discontinuity{$T}:\n   ", d)
