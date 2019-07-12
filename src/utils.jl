"""
    has_constant_lags(integrator::DDEIntegrator)

Return if the DDE problem of the `integrator` contains constant delays.
"""
has_constant_lags(integrator::DDEIntegrator) = has_constant_lags(integrator.sol.prob)

"""
    has_dependent_lags(integrator::DDEIntegrator)

Return if the DDE problem of the `integrator` contains dependent delays.
"""
has_dependent_lags(integrator::DDEIntegrator) = has_dependent_lags(integrator.sol.prob)

"""
    has_constant_lags(prob::DDEProblem)

Return if the DDE problem `prob` contains constant delays.
"""
has_constant_lags(prob::DDEProblem) =
  prob.constant_lags !== nothing && !isempty(prob.constant_lags)

"""
    has_dependent_lags(prob::DDEProblem)

Return if the DDE problem `prob` contains dependent delays.
"""
has_dependent_lags(prob::DDEProblem) =
  prob.dependent_lags !== nothing && !isempty(prob.dependent_lags)

"""
    u_uprev_uprev2(prob, alg; kwargs...)

Return state vectors `u`,`uprev`, and `uprev2` (possibly aliased) for solving the
differential equation problem `prob` with algorithm `alg`.
"""
function u_uprev_uprev2(prob, alg;
                        alias_u0 = false,
                        adaptive = DiffEqBase.isadaptive(alg),
                        allow_extrapolation = alg_extrapolates(alg),
                        calck = false)
  if alias_u0
    u = prob.u0
  else
    u = recursivecopy(prob.u0)
  end

  # Some algorithms do not use `uprev` explicitly. In that case, we can save
  # some memory by aliasing `uprev = u`, e.g. for "2N" low storage methods.
  if OrdinaryDiffEq.uses_uprev(alg, adaptive) || calck
    uprev = recursivecopy(u)
  else
    uprev = u
  end

  if allow_extrapolation
    uprev2 = recursivecopy(u)
  else
    uprev2 = uprev
  end

  u, uprev, uprev2
end

"""
    callback_set_and_cache(prob, callback)

Return set of callbacks and its cache for the differential equation problem `prob` and the
user-provided `callback`.
"""
function callback_set_and_cache(prob, callback)
  callback_set = CallbackSet(callback, prob.callback)

  max_len_cb = DiffEqBase.max_vector_callback_length(callback_set)
  if max_len_cb isa VectorContinuousCallback
    uBottomEltype = recursive_bottom_eltype(prob.u0)
    callback_cache = DiffEqBase.CallbackCache(max_len_cb.len, uBottomEltype, uBottomEltype)
  else
    callback_cache = nothing
  end

  callback_set, callback_cache
end

"""
    rate_prototype_of(prob::DEProblem)

Return prototype of rates for a given differential equation problem.
"""
rate_prototype_of(prob::DiffEqBase.DEProblem) =
  DiffEqBase.@.. prob.u0 * $(inv(oneunit(eltype(prob.tspan))))

"""
    solution_arrays(u, tspan, rate_prototype; kwargs...)

Return arrays of saved time points, states, and rates, initialized with the solution at the
first time point if `save_start = true` (the default).
"""
function solution_arrays(u, tspan, rate_prototype;
                         timeseries_init = typeof(u)[],
                         ts_init = eltype(tspan)[],
                         ks_init = [],
                         save_idxs = nothing,
                         save_start = true)
  # determine types of time and state
  uType = typeof(u)
  tType = eltype(tspan)

  # initialize vector of saved time points
  ts = convert(Vector{tType}, ts_init)

  # initialize vector of saved states
  if save_idxs === nothing
    timeseries = convert(Vector{uType}, timeseries_init)
  else
    u_initial = u[save_idxs]
    timeseries = convert(Vector{typeof(u_initial)}, timeseries_init)
  end

  # initialize vector of saved rates
  if save_idxs === nothing
    ksEltype = Vector{typeof(rate_prototype)}
  else
    ks_prototype = rate_prototype[save_idxs]
    ksEltype = Vector{typeof(ks_prototype)}
  end
  ks = convert(Vector{ksEltype}, ks_init)

  # save solution at initial time point
  if save_start
    copyat_or_push!(ts, 1, first(tspan))
    if save_idxs === nothing
      copyat_or_push!(timeseries, 1, u)
      copyat_or_push!(ks, 1, [rate_prototype])
    else
      u_initial = u[save_idxs]
      copyat_or_push!(timeseries, 1, u_initial, Val{false})
      copyat_or_push!(ks, 1, [ks_prototype])
    end
  end

  ts, timeseries, ks
end

"""
    sizehint!(sol::DESolution, n)

Suggest that solution `sol` reserves capacity for at least `n` elements.
"""
function Base.sizehint!(sol::DESolution, n)
  sizehint!(sol.u, n)
  sizehint!(sol.t, n)
  sizehint!(sol.k, n)

  nothing
end

"""
    sizehint!(sol::DESolution, alg, tspan, tstops, saveat; kwargs...)

Suggest that solution `sol` reserves capacity for a number of elements that
depends on the parameter settings of the numerical solver.
"""
function Base.sizehint!(sol::DESolution, alg, tspan, tstops, saveat;
                        save_everystep = isempty(saveat),
                        adaptive = isadaptive(alg),
                        internalnorm = DiffEqBase.ODE_DEFAULT_NORM,
                        dt = zero(eltype(tspan)))
  # obtain integration time
  t0 = first(tspan)
  integrationtime = last(tspan) - t0

  if !adaptive && save_everystep && !isinf(integrationtime)
    # determine number of steps if known a priori
    if iszero(dt)
      steps = length(tstops)
    else
      steps = ceil(Int, internalnorm(integrationtime / dt, t0))
    end

    sizehint!(sol, steps + 1)
  elseif save_everystep
    sizehint!(sol, 50)
  elseif !isempty(saveat)
    sizehint!(sol, length(saveat) + 1)
  else
    sizehint!(sol, 2)
  end

  nothing
end

"""
    initialize_solution!(integrator)

Initialize the solution of an integrator by adjusting the cache for composite algorithms.
"""
initialize_solution!(integrator) = _initialize_solution!(integrator)
function initialize_solution!(integrator::DDEIntegrator)
  initialize_solution!(integrator.integrator)
  _initialize_solution!(integrator)

  nothing
end

function _initialize_solution!(integrator)
  if integrator.opts.save_start && iscomposite(integrator.alg)
    copyat_or_push!(integrator.sol.alg_choice, 1, integrator.cache.current)
  end

  nothing
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
    :(A = f.f1; expA = exp(A*dt))
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
assign_expr(::Val{name}, ::Type{<:ForwardDiff.JacobianConfig},
            ::Type) where {name} =
                :($name = OrdinaryDiffEq.build_jac_config(alg, f, uf, du1,
                                                            uprev, u, tmp, dz))
assign_expr(::Val{name}, ::Type{<:ForwardDiff.JacobianConfig},
            ::Type{<:OrdinaryDiffEq.RosenbrockMutableCache}) where {name} =
                :($name = OrdinaryDiffEq.build_jac_config(alg, f, uf, du1,
                                                            uprev, u, tmp, du2))
assign_expr(::Val{name}, ::Type{<:DiffEqDiffTools.JacobianCache},
            ::Type) where {name} =
                :($name = OrdinaryDiffEq.build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz))
assign_expr(::Val{name}, ::Type{<:DiffEqDiffTools.JacobianCache},
            ::Type{<:OrdinaryDiffEq.RosenbrockMutableCache}) where {name} =
                :($name = OrdinaryDiffEq.build_jac_config(alg, f,
                            uf, du1, uprev, u, tmp, du2))
# create new config of Gradient
assign_expr(::Val{name}, ::Type{<:ForwardDiff.DerivativeConfig},
           ::Type) where {name} =
               :($name = OrdinaryDiffEq.build_grad_config(alg, f, tf, du1, t))
assign_expr(::Val{name}, ::Type{<:DiffEqDiffTools.GradientCache},
          ::Type) where {name} =
              :($name = OrdinaryDiffEq.build_grad_config(alg, f, tf, du1, t))

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

function unwrap_alg(integrator::DDEIntegrator, is_stiff)
  alg = integrator.alg
  iscomp = typeof(alg) <: CompositeAlgorithm
  if !iscomp
    return alg
  elseif typeof(alg.choice_function) <: AutoSwitch
    num = is_stiff ? 2 : 1
    return alg.algs[num]
  else
    return alg.algs[integrator.cache.current]
  end
end


DiffEqBase.nlsolve_f(integrator::DDEIntegrator) =
  DiffEqBase.nlsolve_f(integrator.f, unwrap_alg(integrator, true))
