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
function compute_discontinuity_tree(lags, alg, start_val,end_val,neutral)
  isempty(lags) && return eltype(lags)[]
  if !neutral
    start_val + unique(vcat((sum.(collect(with_replacement_combinations(lags, i))) for i in
                             1:alg_order(alg))...))
  else
  # Find the minimal amount that the tree must be expanded out to cover the
  # full interval
  # tspan0 + minlag*mult > tspan1 => mult > (tspan1-tspan0)/minlag
    mult = round(Int,(end_val-start_val)/minimum(lags),RoundUp)
    start_val + unique(vcat((sum.(collect(with_replacement_combinations(lags, i))) for i in
                             1:mult)...))
  end
end

"""
    build_linked_cache(cache, alg, u, uprev, uprev2, f, t, dt)

Create cache for algorithm `alg` from existing cache `cache` with updated `u`, `uprev`,
`uprev2`, `f`, `t`, and `dt`.
"""
@generated function build_linked_cache(cache, alg, u, uprev, uprev2, f, t, dt)
    assignments = [assign_expr(Val{name}(), fieldtype(cache, name), cache)
                   for name in fieldnames(cache) if name ∉ [:u, :uprev, :uprev2, :t, :dt]]

    :($(assignments...); $(DiffEqBase.parameterless_type(cache))($(fieldnames(cache)...)))
end

"""
    assign_expr(::Val{name}, ::Type{T}, ::Type{cache})

Create expression that extracts field `name` of type `T` from cache of type `cache`
to variable `name`.

Hereby u, uprev, uprev2, and function f are updated, if required.
"""
assign_expr(::Val{name}, ::Type, ::Type) where {name} =
    :($name = getfield(cache, $(Meta.quot(name))))

# update uhold
assign_expr(::Val{:uhold}, ::Type,
            ::Type{<:Union{OrdinaryDiffEq.GenericImplicitEulerCache,
                           OrdinaryDiffEq.GenericTrapezoidCache,
                           OrdinaryDiffEq.GenericIIF1Cache,
                           OrdinaryDiffEq.GenericIIF2Cache}}) =
                               :(uhold = vec(u))

# update matrix exponential
assign_expr(::Val{:expA}, ::Type, ::Type) =
    :(A = f.f1; expA = expm(A*dt))
assign_expr(::Val{:phi1}, ::Type, ::Type{<:OrdinaryDiffEq.NorsettEulerCache}) =
    :(phi1 = ((expA-I)/A))

# update derivative wrappers
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.TimeDerivativeWrapper}, ::Type) where name =
    :($name = OrdinaryDiffEq.TimeDerivativeWrapper(f, u))
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.UDerivativeWrapper}, ::Type) where name =
    :($name = OrdinaryDiffEq.UDerivativeWrapper(f, getfield(cache, t)))
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.TimeGradientWrapper}, ::Type) where name =
    :($name = OrdinaryDiffEq.TimeGradientWrapper(
        OrdinaryDiffEq.VectorF(f, size(u)),
        uprev,
        getfield(cache, $(Meta.quot(name))).fx1))
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.UJacobianWrapper}, ::Type) where name =
    :($name = OrdinaryDiffEq.UJacobianWrapper(
        OrdinaryDiffEq.VectorFReturn(f, size(u)),
        t,
        vec(uprev),
        getfield(cache, $(Meta.quot(name))).fx1))

# create new config of Jacobian
assign_expr(::Val{name}, ::Type{ForwardDiff.JacobianConfig{T,V,N,D}},
            ::Type) where {name,T,V,N,D} =
                :($name = ForwardDiff.JacobianConfig(uf, vec(du1), vec(uprev),
                                                     ForwardDiff.Chunk{$N}()))

# update implicit RHS
function assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.ImplicitRHS}, ::Type) where name
    nameq = Meta.quot(name)
    :($name = OrdinaryDiffEq.ImplicitRHS(
        f,
        getfield(cache, $nameq).C,
        t, t, t,
        getfield(cache, $nameq).dual_cache))
end
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.ImplicitRHS_Scalar}, ::Type) where name =
    :($name = OrdinaryDiffEq.ImplicitRHS_Scalar(
        f,
        getfield(cache, $(Meta.quot(name))).C,
        t, t, t))
function assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.RHS_IIF}, ::Type) where name
    nameq = Meta.quot(name)
    :($name = OrdinaryDiffEq.RHS_IIF(
        f,
        getfield(cache, $nameq).tmp,
        t, t,
        getfield(cache, $nameq).dual_cache,
        getfield(cache, $nameq).a))
end
function assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.RHS_IIF_Scalar},
                     ::Type) where name
    nameq = Meta.quot(name)
    :($name = OrdinaryDiffEq.RHS_IIF_Scalar(
        f,
        t, t,
        getfield(cache, $nameq).tmp,
        getfield(cache, $nameq).a))
end

# create new NLsolve differentiable function
assign_expr(::Val{name}, ::Type{<:NLsolve.DifferentiableMultivariateFunction},
            ::Type) where name =
                :($name = alg.nlsolve(Val{:init},rhs,uhold))
