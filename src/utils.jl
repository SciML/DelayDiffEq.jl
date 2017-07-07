"""
    compute_discontinuity_tree(lags, alg, start_val)

Compute time points of discontinuities for delays `lags` up to the order of algorithm `alg`,
starting at initial time point `start_val`.

# Examples

```jldoctest
julia> compute_discontinuity_tree([1//2], BS3(), 1)
3-element Array{Rational{Int64},1}:
 3//2
 2//1
 5//2

julia> compute_discontinuity_tree([1//2, 1//3], BS3(), 1)
8-element Array{Rational{Int64},1}:
  3//2
  4//3
  2//1
 11//6
  5//3
  5//2
  7//3
 13//6
```
"""
function compute_discontinuity_tree(lags, alg, start_val)
    start_val + unique(vcat((sum.(collect(with_replacement_combinations(lags, i))) for i in
                             1:alg_order(alg))...))
end
