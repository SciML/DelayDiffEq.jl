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
    u_uprev(u0, alg; kwargs...)

Return state vectors `u` and `uprev` (possibly aliased) for solving the
differential equation problem for initial state `u0` with algorithm `alg`.
"""
function u_uprev(u0, alg;
                 alias_u0 = false,
                 adaptive = DiffEqBase.isadaptive(alg),
                 calck = false)
  if alias_u0
    u = u0
  else
    u = recursivecopy(u0)
  end

  # Some algorithms do not use `uprev` explicitly. In that case, we can save
  # some memory by aliasing `uprev = u`, e.g. for "2N" low storage methods.
  if OrdinaryDiffEq.uses_uprev(alg, adaptive) || calck
    uprev = recursivecopy(u)
  else
    uprev = u
  end

  u, uprev
end

"""
    u_uprev_uprev2(u0, alg; kwargs...)

Return state vectors `u`, `uprev`, and `uprev2` (possibly aliased) for solving the
differential equation problem for initial state `u0` with algorithm `alg`.
"""
function u_uprev_uprev2(u0, alg;
                        allow_extrapolation = alg_extrapolates(alg),
                        kwargs...)
  # compute u and uprev first
  u, uprev = u_uprev(u0, alg; kwargs...)

  if allow_extrapolation
    uprev2 = recursivecopy(u)
  else
    uprev2 = uprev
  end

  u, uprev, uprev2
end

"""
    get_abstol(u, tspan, alg; abstol = nothing)

Return the absolute tolerance for solving the differential equation problem with state
variable `u` and time span `tspan` with algorithm `alg`.
"""
function get_abstol(u, tspan, alg; abstol = nothing)
  if typeof(alg) <: FunctionMap
    _abstol = real.(zero.(u))
  elseif abstol === nothing
    uBottomEltype = recursive_bottom_eltype(u)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

    if uBottomEltypeNoUnits == uBottomEltype
      _abstol = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1//10^6))
    else
      _abstol = real.(oneunit.(u).*1//10^6)
    end
  else
    _abstol = real.(abstol)
  end

  _abstol
end

"""
    get_reltol(u, tspan, alg; reltol = nothing)

Return the relative tolerance for solving the differential equation problem with state
variable `u` and time span `tspan` with algorithm `alg`.
"""
function get_reltol(u, tspan, alg; reltol = nothing)
  if typeof(alg) <: FunctionMap
    _reltol = real.(zero(first(u)/t))
  elseif reltol === nothing
    uBottomEltype = recursive_bottom_eltype(u)
    uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

    if uBottomEltypeNoUnits == uBottomEltype
      _reltol = real(convert(uBottomEltype, oneunit(uBottomEltype) * 1//10^3))
    else
      _reltol = real.(oneunit.(u).*1//10^3)
    end
  else
    _reltol = real.(reltol)
  end

  _reltol
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
    rate_prototype_of(u0, tspan)

Return prototype of rates for a given differential equation problem with state `u` and
time span `tspan`.
"""
rate_prototype_of(u0, tspan) =
  DiffEqBase.@.. u0 * $(inv(oneunit(eltype(tspan))))

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

function build_history_function(prob, alg, rate_prototype, reltol;
                                dt = zero(eltype(prob.tspan)),
                                adaptive = DiffEqBase.isadaptive(alg.alg),
                                calck = false,
                                internalnorm = DiffEqBase.ODE_DEFAULT_NORM)
  @unpack f, u0, tspan, p = prob

  t0 = first(tspan)
  tType = eltype(tspan)
  tTypeNoUnits = typeof(one(tType))
  tdir = sign(last(tspan) - t0)

  uEltypeNoUnits = recursive_unitless_eltype(u0)
  uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u0)

  # bootstrap an ODE integrator
  # - whose solution captures the dense history of the simulation
  # - that is used for extrapolation of the history for time points past the
  #   already fixed history
  # - that is used for interpolation of the history for time points in the
  #   current integration step (so the interpolation is fixed while updating the stages)
  # we wrap the user-provided history function such that function calls during the setup
  # of the integrator do not fail
  ode_f = ODEFunctionWrapper(f, prob.h)
  ode_prob = ODEProblem{isinplace(prob)}(ode_f, u0, tspan, p)

  # get states of ODE integrator (do not alias uprev)
  ode_u, ode_uprev = u_uprev(u0, alg; alias_u0 = false, calck = true)

  # initialize output arrays
  ode_k = typeof(rate_prototype)[]
  ode_ts, ode_timeseries, ode_ks = solution_arrays(ode_u, tspan, rate_prototype;
                                                   save_idxs = nothing,
                                                   save_start = true)

  # obtain cache (we alias uprev2 and uprev)
  ode_cache = OrdinaryDiffEq.alg_cache(alg.alg, ode_u, rate_prototype, uEltypeNoUnits,
                                       uBottomEltypeNoUnits, tTypeNoUnits, ode_uprev,
                                       ode_uprev, ode_f, t0, zero(tType), reltol, p, calck,
                                       Val{isinplace(prob)})

  # build dense interpolation of history
  if iscomposite(alg)
    ode_alg_choice = Int[]
    ode_id = OrdinaryDiffEq.CompositeInterpolationData(ode_f, ode_timeseries, ode_ts, ode_ks,
                                                       ode_alg_choice, true, ode_cache) # dense = true
    ode_sol = DiffEqBase.build_solution(ode_prob, alg.alg, ode_ts, ode_timeseries;
                                        dense = true, k = ode_ks, interp = ode_id,
                                        alg_choice = ode_alg_choice,
                                        calculate_error = false, destats = DiffEqBase.DEStats(0))
  else
    ode_id = OrdinaryDiffEq.InterpolationData(ode_f, ode_timeseries, ode_ts, ode_ks, true, ode_cache) # dense = true
    ode_sol = DiffEqBase.build_solution(ode_prob, alg.alg, ode_ts, ode_timeseries;
                                        dense = true, k = ode_ks, interp = ode_id,
                                        calculate_error = false, destats = DiffEqBase.DEStats(0))
  end

  # reserve capacity
  sizehint!(ode_sol, alg.alg, tspan, (), ();
            save_everystep = true, adaptive = adaptive, internalnorm = internalnorm, dt = tType(dt))

  # create simple integrator
  tdirType = typeof(sign(zero(tType)))
  ode_integrator = HistoryODEIntegrator{typeof(alg.alg),isinplace(prob),typeof(prob.u0),
                                        tType,tdirType,typeof(ode_k),
                                        typeof(ode_sol),typeof(ode_cache)}(
                                          ode_sol, ode_u, ode_k, t0, zero(tType), ode_uprev,
                                          t0, alg.alg, zero(tType), tdir, 1, 1, ode_cache)

  # combine the user-provided history function and the ODE integrator with dense solution
  # to a joint dense history of the DDE
  # we use this history information to create a problem function of the DDE with all
  # available history information that is of the form f(du,u,p,t) or f(u,p,t) such that
  # ODE algorithms can be applied
  HistoryFunction(prob.h, ode_integrator)
end

"""
    initialize_solution!(integrator::DDEIntegrator)

Initialize the solution of an integrator by adjusting the cache for composite algorithms.
"""
function initialize_solution!(integrator::DDEIntegrator)
  if iscomposite(integrator.alg)
    copyat_or_push!(integrator.integrator.sol.alg_choice, 1, integrator.cache.current)
    if integrator.opts.save_start
      copyat_or_push!(integrator.sol.alg_choice, 1, integrator.cache.current)
    end
  end

  nothing
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
