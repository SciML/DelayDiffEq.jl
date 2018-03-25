"""
    agrees(h, u, p, t)

Determine whether history function evaluates to `u` at time point `t`
for parameters `p`.
"""
function agrees(h, u, p, t)
    # Obtain signatures of h
    sigs = [m.sig for m in methods(h)]

    # Compare evaluation of h at time point t with u
    if any(sig<:Tuple{Any, Any, Any} for sig in sigs)
        return h(p, t) == u
    elseif any(sig<:Tuple{Any, Any, Any, Any} for sig in sigs)
        val = recursivecopy(u)
        h(val, p, t)
        return val == u
    end

    return false
end


"""
    fsal_typeof(integrator::ODEIntegrator)

Return type of FSAL of `integrator`.
"""
function fsal_typeof(integrator::ODEIntegrator{<:OrdinaryDiffEq.OrdinaryDiffEqAlgorithm,
                                               uType,tType,P,tTypeNoUnits,tdirType,ksEltype,
                                               SolType,F,ProgressType,CacheType,O,
                                               FSALType}) where {uType,tType,P,tTypeNoUnits,
                                                                 tdirType,ksEltype,SolType,
                                                                 F,ProgressType,CacheType,O,
                                                                 FSALType}
    return FSALType
end

"""
    build_linked_cache(cache, alg, u, uprev, uprev2, f, t, dt)

Create cache for algorithm `alg` from existing cache `cache` with updated `u`, `uprev`,
`uprev2`, `f`, `t`, and `dt`.
"""
@generated function build_linked_cache(cache, alg, u, uprev, uprev2, f, t, dt,p)
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

# update matrix exponential
assign_expr(::Val{:expA}, ::Type, ::Type) =
    :(A = f.f1; expA = expm(A*dt))
assign_expr(::Val{:phi1}, ::Type, ::Type{<:OrdinaryDiffEq.NorsettEulerCache}) =
    :(phi1 = ((expA-I)/A))

# update derivative wrappers
assign_expr(::Val{name}, ::Type{<:DiffEqDiffTools.TimeDerivativeWrapper}, ::Type) where name =
    :($name = DiffEqDiffTools.TimeDerivativeWrapper(f, u,p))
assign_expr(::Val{name}, ::Type{<:DiffEqDiffTools.UDerivativeWrapper}, ::Type) where name =
    :($name = DiffEqDiffTools.UDerivativeWrapper(f, t,p))
assign_expr(::Val{name}, ::Type{<:DiffEqDiffTools.TimeGradientWrapper}, ::Type) where name =
    :($name = DiffEqDiffTools.TimeGradientWrapper(
        f,uprev,p))
assign_expr(::Val{name}, ::Type{<:DiffEqDiffTools.UJacobianWrapper}, ::Type) where name =
    :($name = DiffEqDiffTools.UJacobianWrapper(
        f,t,p))

# create new config of Jacobian
assign_expr(::Val{name}, ::Type{ForwardDiff.JacobianConfig{T,V,N,D}},
            ::Type) where {name,T,V,N,D} =
                :($name = ForwardDiff.JacobianConfig(uf, du1, uprev,
                                                     ForwardDiff.Chunk{$N}()))

# update implicit RHS
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.ImplicitRHS}, ::Type) where name =
    :($name = OrdinaryDiffEq.ImplicitRHS(f, cache.tmp, t, t, t, cache.dual_cache,p))
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.ImplicitRHS_Scalar}, ::Type) where name =
    :($name = OrdinaryDiffEq.ImplicitRHS_Scalar(f, zero(u), t, t, t,p))
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.RHS_IIF}, ::Type) where name =
    :($name = OrdinaryDiffEq.RHS_IIF(f, cache.tmp, t, t, cache.tmp, cache.dual_cache,p))
assign_expr(::Val{name}, ::Type{<:OrdinaryDiffEq.RHS_IIF_Scalar}, ::Type) where name =
    :($name = OrdinaryDiffEq.RHS_IIF_Scalar(f, zero(u), t, t,
                                            getfield(cache, $(Meta.quot(name))).a,p))

# create new NLsolve differentiable function
assign_expr(::Val{name}, ::Type{<:NLSolversBase.OnceDifferentiable},
            ::Type{<:OrdinaryDiffEq.OrdinaryDiffEqMutableCache}) where name =
                :($name = alg.nlsolve(Val{:init},rhs,u))
assign_expr(::Val{name}, ::Type{<:NLSolversBase.OnceDifferentiable},
            ::Type{<:OrdinaryDiffEq.OrdinaryDiffEqConstantCache}) where name =
                :($name = alg.nlsolve(Val{:init},rhs,uhold))
