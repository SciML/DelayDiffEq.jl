function DiffEqBase.__solve(prob::DiffEqBase.AbstractDDEProblem,
                            alg::AbstractMethodOfStepsAlgorithm, args...;
                            kwargs...)
  integrator = DiffEqBase.__init(prob, alg, args...; kwargs...)
  DiffEqBase.solve!(integrator)
  integrator.sol
end

function DiffEqBase.__init(prob::DiffEqBase.AbstractDDEProblem,
                           alg::AbstractMethodOfStepsAlgorithm,
                           timeseries_init = typeof(prob.u0)[],
                           ts_init = eltype(prob.tspan)[],
                           ks_init = [];
                           saveat = eltype(prob.tspan)[],
                           tstops = eltype(prob.tspan)[],
                           d_discontinuities = Discontinuity{eltype(prob.tspan)}[],
                           save_idxs = nothing,
                           save_everystep = isempty(saveat),
                           save_on = true,
                           save_start = save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[1] in saveat,
                           save_end = save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[end] in saveat,
                           callback = nothing,
                           dense = save_everystep && isempty(saveat),
                           calck = (callback !== nothing && callback != CallbackSet()) || # Empty callback
                                   (!isempty(setdiff(saveat,tstops)) || dense), # and no dense output
                           dt = zero(eltype(prob.tspan)),
                           dtmin = typeof(one(eltype(prob.tspan))) <: AbstractFloat ? eps(eltype(prob.tspan)) :
                                   typeof(one(eltype(prob.tspan))) <: Integer ? 0 : eltype(prob.tspan)(1//10^(10)),
                           dtmax = eltype(prob.tspan)(prob.tspan[end]-prob.tspan[1]),
                           force_dtmin = false,
                           adaptive = DiffEqBase.isadaptive(alg),
                           gamma = OrdinaryDiffEq.gamma_default(alg.alg),
                           abstol = nothing,
                           reltol = nothing,
                           qmin = OrdinaryDiffEq.qmin_default(alg.alg),
                           qmax = OrdinaryDiffEq.qmax_default(alg.alg),
                           qsteady_min = OrdinaryDiffEq.qsteady_min_default(alg.alg),
                           qsteady_max = OrdinaryDiffEq.qsteady_max_default(alg.alg),
                           qoldinit = 1//10^4,
                           fullnormalize = true,
                           failfactor = 2,
                           beta1 = nothing,
                           beta2 = nothing,
                           maxiters = adaptive ? 1000000 : typemax(Int),
                           internalnorm = DiffEqBase.ODE_DEFAULT_NORM,
                           internalopnorm = LinearAlgebra.opnorm,
                           isoutofdomain = DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN,
                           unstable_check = DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK,
                           verbose = true,
                           timeseries_errors = true,
                           dense_errors = false,
                           advance_to_tstop = false,
                           stop_at_next_tstop = false,
                           initialize_save = true,
                           progress = false,
                           progress_steps = 1000,
                           progress_name = "DDE",
                           progress_message = DiffEqBase.ODE_DEFAULT_PROG_MESSAGE,
                           userdata = nothing,
                           allow_extrapolation = OrdinaryDiffEq.alg_extrapolates(alg),
                           initialize_integrator = true,
                           alias_u0 = false,
                           # keyword arguments for DDEs
                           discontinuity_interp_points::Int = 10,
                           discontinuity_abstol = eltype(prob.tspan)(1//Int64(10)^12),
                           discontinuity_reltol = 0,
                           kwargs... )
  if haskey(kwargs, :initial_order)
    @warn "initial_order has been deprecated. Please specify order_discontinuity_t0 in the DDEProblem instead."
    order_discontinuity_t0 = kwargs[:initial_order]
  else
    order_discontinuity_t0 = prob.order_discontinuity_t0
  end

  if haskey(kwargs, :minimal_solution)
    @warn "minimal_solution is ignored"
  end

  if !isempty(saveat) && dense
    @warn("Dense output is incompatible with saveat. Please use the SavingCallback from the Callback Library to mix the two behaviors.")
  end

  progress && @logmsg(-1, progress_name, _id = _id = :OrdinaryDiffEq, progress = 0)

  # unpack problem
  @unpack f, u0, h, tspan, p, neutral, constant_lags, dependent_lags = prob

  # determine type and direction of time
  tType = eltype(tspan)
  t0 = first(tspan)
  tdir = sign(last(tspan) - t0)
  tTypeNoUnits = typeof(one(tType))

  # Allow positive dtmax, but auto-convert
  dtmax > zero(dtmax) && tdir < zero(tdir) && (dtmax *= tdir)

  # no fixed-point iterations for constrained algorithms,
  # and thus `dtmax` should match minimal lag
  if isconstrained(alg) && has_constant_lags(prob)
    dtmax = tdir * min(abs(dtmax), minimum(abs, constant_lags))
  end

  # get absolute and relative tolerances
  abstol_internal = get_abstol(u0, tspan, alg.alg; abstol = abstol)
  reltol_internal = get_reltol(u0, tspan, alg.alg; reltol = reltol)

  # get rate prototype
  rate_prototype = rate_prototype_of(u0, tspan)

  # create a history function
  history = build_history_function(prob, alg, rate_prototype, reltol_internal;
                                   dt = dt, adaptive = adaptive,
                                   internalnorm = internalnorm)
  f_with_history = ODEFunctionWrapper(f, history)

  # get states (possibly different from the ODE integrator!)
  u, uprev, uprev2 = u_uprev_uprev2(u0, alg;
                                    alias_u0 = alias_u0,
                                    adaptive = adaptive,
                                    allow_extrapolation = allow_extrapolation,
                                    calck = calck)
  uEltypeNoUnits = recursive_unitless_eltype(u)
  uBottomEltypeNoUnits = recursive_unitless_bottom_eltype(u)

  # initialize output arrays of the solution
  k = typeof(rate_prototype)[]
  ts, timeseries, ks = solution_arrays(u, tspan, rate_prototype;
                                       timeseries_init = timeseries_init,
                                       ts_init = ts_init,
                                       ks_init = ks_init,
                                       save_idxs = save_idxs,
                                       save_start = save_start)

  # build cache
  ode_integrator = history.integrator
  cache = OrdinaryDiffEq.alg_cache(alg.alg, u, rate_prototype, uEltypeNoUnits,
                                   uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2,
                                   f_with_history, t0, zero(tType), reltol_internal, p, calck,
                                   Val(isinplace(prob)))

  # separate statistics of the integrator and the history
  destats = DDEStats(0)

  # create solution
  if iscomposite(alg)
    id = OrdinaryDiffEq.CompositeInterpolationData(f_with_history, timeseries, ts, ks,
                                                   Int[], dense, cache)
    sol = DiffEqBase.build_solution(prob, alg.alg, ts, timeseries;
                                    dense = dense, k = ks, interp = id,
                                    alg_choice = id.alg_choice, calculate_error = false,
                                    destats = destats)
  else
    id = OrdinaryDiffEq.InterpolationData(f_with_history, timeseries, ts, ks, dense, cache)
    sol = DiffEqBase.build_solution(prob, alg.alg, ts, timeseries;
                                    dense = dense, k = ks, interp = id,
                                    calculate_error = false, destats = destats)
  end

  # retrieve time stops, time points at which solutions is saved, and discontinuities
  maximum_order = OrdinaryDiffEq.alg_maximum_order(alg)
  tstops_internal, saveat_internal, d_discontinuities_internal =
    OrdinaryDiffEq.tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan,
                                              order_discontinuity_t0, maximum_order,
                                              constant_lags, neutral)

  # reserve capacity for the solution
  sizehint!(sol, alg, tspan, tstops_internal, saveat_internal;
            save_everystep = save_everystep, adaptive = adaptive, dt = tType(dt))

  # create array of tracked discontinuities
  # used to find propagated discontinuities with callbacks and to keep track of all
  # passed discontinuities
  tracked_discontinuities = Discontinuity{tType}[]
  if order_discontinuity_t0 ≤ maximum_order
    push!(tracked_discontinuities, Discontinuity(tdir * t0, order_discontinuity_t0))
  end

  # Create set of callbacks and its cache
  callback_set, callback_cache = callback_set_and_cache(prob, callback)

  # separate options of integrator and of dummy ODE integrator since ODE integrator always saves
  # every step and every index (necessary for history function)
  QT = tTypeNoUnits <: Integer ? typeof(qmin) : tTypeNoUnits

  if iscomposite(alg)
    if beta2 === nothing
      _beta2 = QT(OrdinaryDiffEq.beta2_default(alg.alg.algs[cache.current]))
    else
      _beta2 = QT(beta2)
    end

    if beta1 === nothing
      _beta1 = QT(OrdinaryDiffEq.beta1_default(alg.alg.algs[cache.current], _beta2))
    else
      _beta1 = QT(beta1)
    end
  else
    if beta2 === nothing
      _beta2 = QT(OrdinaryDiffEq.beta2_default(alg.alg))
    else
      _beta2 = QT(beta2)
    end

    if beta1 === nothing
      _beta1 = QT(OrdinaryDiffEq.beta1_default(alg.alg, _beta2))
    else
      _beta1 = QT(beta1)
    end
  end

  opts = OrdinaryDiffEq.DEOptions{typeof(abstol_internal),typeof(reltol_internal),QT,tType,
                                  typeof(internalnorm),typeof(internalopnorm),
                                  typeof(callback_set),typeof(isoutofdomain),
                                  typeof(progress_message),typeof(unstable_check),
                                  typeof(tstops_internal),
                                  typeof(d_discontinuities_internal),typeof(userdata),
                                  typeof(save_idxs),typeof(maxiters),typeof(tstops),
                                  typeof(saveat),typeof(d_discontinuities)}(
                                    maxiters,save_everystep,adaptive,abstol_internal,
                                    reltol_internal,QT(gamma),QT(qmax),
                                    QT(qmin),QT(qsteady_max),
                                    QT(qsteady_min),QT(failfactor),tType(dtmax),
                                    tType(dtmin),internalnorm,internalopnorm,save_idxs,
                                    tstops_internal,saveat_internal,
                                    d_discontinuities_internal,
                                    tstops,saveat,d_discontinuities,
                                    userdata,progress,progress_steps,
                                    progress_name,progress_message,timeseries_errors,
                                    dense_errors,_beta1,_beta2,QT(qoldinit),dense,
                                    save_on,save_start,save_end,callback_set,
                                    isoutofdomain,unstable_check,verbose,calck,force_dtmin,
                                    advance_to_tstop,stop_at_next_tstop)

  # create fixed point solver
  fpsolver = build_fpsolver(alg, alg.fpsolve, u, uEltypeNoUnits, uBottomEltypeNoUnits,
                            Val(isinplace(prob)))

  # initialize indices of u(t) and u(tprev) in the dense history
  prev_idx = 1
  prev2_idx = 1

  # create integrator combining the new defined problem function with history
  # information, the new solution, the parameters of the ODE integrator, and
  # parameters of fixed-point iteration
  # do not initialize fsalfirst and fsallast
  # rate/state = (state/time)/state = 1/t units, internalnorm drops units
  eigen_est = one(uBottomEltypeNoUnits)/one(tType)
  tprev = t0
  dtcache = tType(dt)
  dtpropose = tType(dt)
  iter = 0
  kshortsize = 0
  reeval_fsal = false
  u_modified = false
  EEst = tTypeNoUnits(1)
  just_hit_tstop = false
  isout = false
  accept_step = false
  force_stepfail = false
  last_stepfail = false
  event_last_time = 0
  vector_event_last_time = 1
  last_event_error = zero(uBottomEltypeNoUnits)
  dtchangeable = OrdinaryDiffEq.isdtchangeable(alg.alg)
  q11 = tTypeNoUnits(1)
  success_iter = 0
  erracc = tTypeNoUnits(1)
  dtacc = tType(1)

  integrator = DDEIntegrator{typeof(alg.alg),isinplace(prob),typeof(u0),tType,typeof(p),
                             typeof(eigen_est),QT,typeof(tdir),typeof(k),typeof(sol),
                             typeof(f_with_history),typeof(cache),
                             typeof(ode_integrator),typeof(fpsolver),
                             typeof(opts),typeof(discontinuity_abstol),
                             typeof(discontinuity_reltol),typeof(history),
                             OrdinaryDiffEq.fsal_typeof(alg.alg, rate_prototype),
                             typeof(last_event_error),typeof(callback_cache)}(
                               sol, u, k, t0, tType(dt), f_with_history, p,
                               uprev, uprev2, tprev, prev_idx, prev2_idx,
                               fpsolver,
                               order_discontinuity_t0, tracked_discontinuities,
                               discontinuity_interp_points, discontinuity_abstol,
                               discontinuity_reltol, alg.alg,
                               dtcache, dtchangeable, dtpropose, tdir,
                               eigen_est, EEst, QT(qoldinit), q11,
                               erracc, dtacc, success_iter,
                               iter, length(ts), length(ts), cache, callback_cache,
                               kshortsize, force_stepfail, last_stepfail,
                               just_hit_tstop, event_last_time, vector_event_last_time,
                               last_event_error, accept_step,
                               isout, reeval_fsal,
                               u_modified, opts, destats, history, ode_integrator)

  # initialize DDE integrator
  if initialize_integrator
    initialize_solution!(integrator)
    OrdinaryDiffEq.initialize_callbacks!(integrator, initialize_save)
    OrdinaryDiffEq.initialize!(integrator)
  end

  # take care of time step dt = 0 and dt with incorrect sign
  OrdinaryDiffEq.handle_dt!(integrator)

  integrator
end

function DiffEqBase.solve!(integrator::DDEIntegrator)
  @unpack tdir, opts, sol = integrator
  @unpack tstops = opts

  # step over all stopping time points, similar to solving with ODE integrators
  @inbounds while !isempty(tstops)
    while tdir * integrator.t < top(tstops)
      # apply step or adapt step size
      OrdinaryDiffEq.loopheader!(integrator)

      # abort integration following same criteria as for ODEs:
      # maxiters exceeded, dt <= dtmin, integration unstable
      DiffEqBase.check_error!(integrator) === :Success || return sol

      # calculate next step
      OrdinaryDiffEq.perform_step!(integrator)

      # calculate proposed next step size, handle callbacks, and update solution
      OrdinaryDiffEq.loopfooter!(integrator)

      isempty(tstops) && break
    end

    # remove hit or passed stopping time points
    OrdinaryDiffEq.handle_tstop!(integrator)
  end

  # clean up solution
  DiffEqBase.postamble!(integrator)

  f = sol.prob.f

  if DiffEqBase.has_analytic(f)
    DiffEqBase.calculate_solution_errors!(sol;
                                          timeseries_errors = opts.timeseries_errors,
                                          dense_errors = opts.dense_errors)
  end
  sol.retcode === :Default || return sol

  integrator.sol = DiffEqBase.solution_new_retcode(sol, :Success)
end

function OrdinaryDiffEq.initialize_callbacks!(integrator::DDEIntegrator,
                                              initialize_save = true)
  callbacks = integrator.opts.callback
  prob = integrator.sol.prob

  # set up additional initial values of newly created DDE integrator
  # (such as fsalfirst) and its callbacks

  integrator.u_modified = true

  u_modified = initialize!(callbacks, integrator.u, integrator.t, integrator)

  # if the user modifies u, we need to fix previous values before initializing
  # FSAL in order for the starting derivatives to be correct
  if u_modified

    if isinplace(prob)
      recursivecopy!(integrator.uprev, integrator.u)
    else
      integrator.uprev = integrator.u
    end

    if OrdinaryDiffEq.alg_extrapolates(integrator.alg)
      if isinplace(prob)
        recursivecopy!(integrator.uprev2, integrator.uprev)
      else
        integrator.uprev2 = integrator.uprev
      end
    end

    # update heap of discontinuities
    # discontinuity is assumed to be of order 0, i.e. solution x is discontinuous
    push!(integrator.opts.d_discontinuities, Discontinuity(integrator.tdir * integrator.t, 0))

    # reset this as it is now handled so the integrators should proceed as normal
    reeval_internals_due_to_modification!(integrator, Val{false})

    if initialize_save &&
      (any((c)->c.save_positions[2],callbacks.discrete_callbacks) ||
       any((c)->c.save_positions[2],callbacks.continuous_callbacks))
      savevalues!(integrator, true)
    end
  end

  # reset this as it is now handled so the integrators should proceed as normal
  integrator.u_modified = false
end

function OrdinaryDiffEq.tstop_saveat_disc_handling(tstops, saveat, d_discontinuities, tspan,
                                                   order_discontinuity_t0, alg_maximum_order, constant_lags,
                                                   neutral)
  tType = eltype(tspan)
  t0, tf = tspan
  tdir = sign(tf - t0)
  tdir_t0 = tdir * t0
  tdir_tf = tdir * tf

  # add discontinuities propagated from initial discontinuity
  add_propagated_constant_lags = order_discontinuity_t0 ≤ alg_maximum_order &&
    constant_lags !== nothing && !isempty(constant_lags)

  if add_propagated_constant_lags
    maxlag = tdir_tf - tdir_t0
    next_order = neutral ? order_discontinuity_t0 : order_discontinuity_t0 + 1
  end

  # time stops
  tstops_internal = BinaryMinHeap{tType}()
  if isempty(d_discontinuities) && !add_propagated_constant_lags && isempty(tstops) # TODO: Specialize more
    push!(tstops_internal, tdir_tf)
  else
    for t in tstops
      tdir_t = tdir * t
      tdir_t0 < tdir_t ≤ tdir_tf && push!(tstops_internal, tdir_t)
    end

    for t in d_discontinuities
      tdir_t = tdir * t
      tdir_t0 < tdir_t ≤ tdir_tf && push!(tstops_internal, tdir_t)
    end

    # add propagated discontinuities
    if add_propagated_constant_lags
      for lag in constant_lags
        if tdir * lag < maxlag
          push!(tstops_internal, tdir * (t0 + lag))
        end
      end
    end

    push!(tstops_internal, tdir_tf)
  end

  # saving time points
  saveat_internal = BinaryMinHeap{tType}()
  if typeof(saveat) <: Number
    if (t0:saveat:tf)[end] == tf
      for t in (t0 + saveat):saveat:tf
        push!(saveat_internal, tdir * t)
      end
    else
      for t in (t0 + saveat):saveat:(tf - saveat)
        push!(saveat_internal, tdir * t)
      end
    end
  elseif !isempty(saveat)
    for t in saveat
      tdir_t = tdir * t
      tdir_t0 < tdir_t < tdir_tf && push!(saveat_internal, tdir_t)
    end
  end

  # discontinuities
  d_discontinuities_internal = BinaryMinHeap{Discontinuity{tType}}()
  if add_propagated_constant_lags
    sizehint!(d_discontinuities_internal.valtree, length(d_discontinuities) + length(constant_lags))
  else
    sizehint!(d_discontinuities_internal.valtree, length(d_discontinuities))
  end

  for d in d_discontinuities
    tdir_t = tdir * d.t

    if tdir_t0 < tdir_t < tdir_tf && d.order ≤ alg_maximum_order + 1
      push!(d_discontinuities_internal, Discontinuity{tType}(tdir_t, d.order))
    end
  end

  if add_propagated_constant_lags
    for lag in constant_lags
      if tdir * lag < maxlag
        push!(d_discontinuities_internal, Discontinuity{tType}(tdir * (t0 + lag), next_order))
      end
    end
  end

  tstops_internal, saveat_internal, d_discontinuities_internal
end
