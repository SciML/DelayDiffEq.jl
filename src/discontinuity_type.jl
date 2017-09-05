"""
    Discontinuity(t, order::Int)

Object of discontinuity of order `order` at time `t`.
"""
struct Discontinuity{tType}
    t::tType
    order::Int
end

# comparisons enable creation of heaps with objects of type `Discontinuity`
# that are automatically ordered
Base.:<(d1::Discontinuity, d2::Discontinuity) = d1.t < d2.t
Base.:>(d1::Discontinuity, d2::Discontinuity) = d1.t > d2.t

# sorting
Base.isless(d1::Discontinuity, d2::Discontinuity) = isless(d1.t, d2.t)

# conversion
Base.convert(::Type{T}, d::Discontinuity) where {T<:Number} = T(d.t)
