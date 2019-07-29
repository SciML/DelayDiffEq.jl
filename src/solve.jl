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
                                   (prob.callback !== nothing && prob.callback != CallbackSet()) || # Empty prob.callback
                                   (!isempty(setdiff(saveat,tstops)) || dense), # and no dense output
                           dt = zero(eltype(prob.tspan)),
                           dtmax = eltype(prob.tspan)(prob.tspan[end]-prob.tspan[1]),
                           adaptive = DiffEqBase.isadaptive(alg),
                           timeseries_errors = true,
                           dense_errors = false,
                           initialize_save = true,
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

  # unpack problem
  @unpack f, u0, h, tspan, p, neutral, constant_lags, dependent_lags = prob

  # determine type and direction of time
  tType = eltype(tspan)
  tdir = sign(last(tspan) - first(tspan))

  # Allow positive dtmax, but auto-convert
  dtmax > zero(dtmax) && tdir < zero(tdir) && (dtmax *= tdir)

  # no fixed-point iterations for constrained algorithms,
  # and thus `dtmax` should match minimal lag
  if isconstrained(alg) && has_constant_lags(prob)
    dtmax = tdir * min(abs(dtmax), minimum(abs, constant_lags))
  end

  # bootstrap an ODE integrator
  # - whose solution captures the dense history of the simulation
  # - that is used for extrapolation of the history for time points past the
  #   already fixed history
  # - that is used for interpolation of the history for time points in the
  #   current integration step (so the interpolation is fixed while updating the stages)
  # we wrap the user-provided history function such that function calls during the setup
  # of the integrator do not fail
  ode_f = ODEFunctionWrapper(f, h)
  ode_prob = ODEProblem{isinplace(prob)}(ode_f, u0, tspan, p)
  ode_integrator = init(ode_prob, alg.alg; initialize_integrator = false, alias_u0 = false,
                        dt = oneunit(tType), dtmax = dtmax, adaptive = adaptive,
                        dense = true, save_everystep = true, save_start = true,
                        save_end = true, kwargs...)

  # ensure that ODE integrator satisfies tprev + dt == t
  ode_integrator.dt = zero(ode_integrator.dt)
  ode_integrator.dtcache = zero(ode_integrator.dt)

  # combine the user-provided history function, the dense solution of the ODE integrator,
  # and the inter- and extrapolation of the integrator to a joint dense history of the
  # DDE
  # we use this history information to create a problem function of the DDE with all
  # available history information that is of the form f(du,u,p,t) or f(u,p,t) such that
  # ODE algorithms can be applied
  history = HistoryFunction(h, ode_integrator)
  f_with_history = ODEFunctionWrapper(f, history)

  # get states (possibly different from the ODE integrator!)
  u, uprev, uprev2 = u_uprev_uprev2(prob, alg;
                                    alias_u0 = alias_u0,
                                    adaptive = adaptive,
                                    allow_extrapolation = allow_extrapolation,
                                    calck = calck)

  # initialize output arrays of the solution
  rate_prototype = rate_prototype_of(prob)
  k = typeof(rate_prototype)[]
  ts, timeseries, ks = solution_arrays(u, tspan, rate_prototype;
                                       timeseries_init = timeseries_init,
                                       ts_init = ts_init,
                                       ks_init = ks_init,
                                       save_idxs = save_idxs,
                                       save_start = save_start)

  # derive cache for states and function with wrapped dense history from the
  # cache of the ODE integrator
  if iscomposite(alg)
    caches = map((x, y) -> build_linked_cache(x, y, u, uprev, uprev2, f_with_history,
                                              tspan[1], dt, p),
                 ode_integrator.cache.caches, alg.alg.algs)
    cache = OrdinaryDiffEq.CompositeCache(caches, alg.alg.choice_function, 1)
  else
    cache = build_linked_cache(ode_integrator.cache, alg.alg, u, uprev, uprev2,
                               f_with_history, tspan[1], dt, p)
  end

  # separate statistics of the integrator and the history
  destats = DiffEqBase.DEStats(0)

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
            save_everystep = save_everystep, adaptive = adaptive, dt = dt)

  # create array of tracked discontinuities
  # used to find propagated discontinuities with callbacks and to keep track of all
  # passed discontinuities
  if order_discontinuity_t0 ≤ maximum_order
    tracked_discontinuities = [Discontinuity(tspan[1], order_discontinuity_t0)]
  else
    tracked_discontinuities = Discontinuity{tType}[]
  end

  # Create set of callbacks and its cache
  callback_set, callback_cache = callback_set_and_cache(prob, callback)

  # separate options of integrator and of dummy ODE integrator since ODE integrator always saves
  # every step and every index (necessary for history function)
  opts = OrdinaryDiffEq.DEOptions(ode_integrator.opts.maxiters,
                                  save_everystep, adaptive,
                                  ode_integrator.opts.abstol,
                                  ode_integrator.opts.reltol,
                                  ode_integrator.opts.gamma,
                                  ode_integrator.opts.qmax,
                                  ode_integrator.opts.qmin,
                                  ode_integrator.opts.qsteady_max,
                                  ode_integrator.opts.qsteady_min,
                                  ode_integrator.opts.failfactor,
                                  ode_integrator.opts.dtmax,
                                  ode_integrator.opts.dtmin,
                                  ode_integrator.opts.internalnorm,
                                  ode_integrator.opts.internalopnorm,
                                  save_idxs, tstops_internal, saveat_internal,
                                  d_discontinuities_internal, tstops, saveat,
                                  d_discontinuities,
                                  ode_integrator.opts.userdata,
                                  ode_integrator.opts.progress,
                                  ode_integrator.opts.progress_steps,
                                  ode_integrator.opts.progress_name,
                                  ode_integrator.opts.progress_message,
                                  timeseries_errors, dense_errors,
                                  ode_integrator.opts.beta1,
                                  ode_integrator.opts.beta2,
                                  ode_integrator.opts.qoldinit,
                                  dense, save_on, save_start, save_end, callback_set,
                                  ode_integrator.opts.isoutofdomain,
                                  ode_integrator.opts.unstable_check,
                                  ode_integrator.opts.verbose, calck,
                                  ode_integrator.opts.force_dtmin,
                                  ode_integrator.opts.advance_to_tstop,
                                  ode_integrator.opts.stop_at_next_tstop)

  # create container for residuals (has to be unitless)
  uEltypeNoUnits = recursive_unitless_eltype(u)
  if typeof(ode_integrator.u) <: AbstractArray
    resid = similar(ode_integrator.u, uEltypeNoUnits)
  else
    resid = one(uEltypeNoUnits)
  end

  # define absolute tolerance for fixed-point iterations
  if alg.fixedpoint_abstol === nothing
    fixedpoint_abstol_internal = recursivecopy(ode_integrator.opts.abstol)
  else
    fixedpoint_abstol_internal = real.(alg.fixedpoint_abstol)
  end

  # use norm of the ODE integrator if no norm for fixed-point iterations is specified
  if alg.fixedpoint_norm === nothing
    fixedpoint_norm = ode_integrator.opts.internalnorm
  end

  # define relative tolerance for fixed-point iterations
  if alg.fixedpoint_reltol === nothing
    fixedpoint_reltol_internal = recursivecopy(ode_integrator.opts.reltol)
  else
    fixedpoint_reltol_internal = real.(alg.fixedpoint_reltol)
  end

  # initialize indices of u(t) and u(tprev) in the dense history
  prev_idx = 1
  prev2_idx = 1

  # create integrator combining the new defined problem function with history
  # information, the new solution, the parameters of the ODE integrator, and
  # parameters of fixed-point iteration
  # do not initialize fsalfirst and fsallast
  tTypeNoUnits = typeof(one(tType))
  integrator = DDEIntegrator{typeof(alg.alg),isinplace(prob),typeof(u),tType,typeof(p),
                             typeof(ode_integrator.eigen_est),
                             typeof(fixedpoint_abstol_internal),
                             typeof(fixedpoint_reltol_internal),typeof(resid),tTypeNoUnits,
                             typeof(tdir),typeof(k),typeof(sol),typeof(f_with_history),
                             typeof(cache),typeof(ode_integrator),typeof(fixedpoint_norm),
                             typeof(opts),typeof(discontinuity_abstol),
                             typeof(discontinuity_reltol),typeof(history),
                             OrdinaryDiffEq.fsal_typeof(alg.alg, rate_prototype),
                             typeof(ode_integrator.last_event_error),
                             typeof(callback_cache)}(
                               sol, u, k, ode_integrator.t, tType(dt), f_with_history, p,
                               uprev, uprev2, ode_integrator.tprev, prev_idx, prev2_idx,
                               fixedpoint_abstol_internal, fixedpoint_reltol_internal,
                               resid, fixedpoint_norm, alg.max_fixedpoint_iters,
                               order_discontinuity_t0, tracked_discontinuities,
                               discontinuity_interp_points, discontinuity_abstol,
                               discontinuity_reltol, alg.alg,
                               tType(dt), ode_integrator.dtchangeable, tType(dt), tdir,
                               ode_integrator.eigen_est, ode_integrator.EEst,
                               ode_integrator.qold, ode_integrator.q11,
                               ode_integrator.erracc, ode_integrator.dtacc,
                               ode_integrator.success_iter, ode_integrator.iter,
                               length(ts), length(ts), cache, callback_cache,
                               ode_integrator.kshortsize, ode_integrator.force_stepfail,
                               ode_integrator.last_stepfail, ode_integrator.just_hit_tstop,
                               ode_integrator.event_last_time,
                               ode_integrator.vector_event_last_time,
                               ode_integrator.last_event_error, ode_integrator.accept_step,
                               ode_integrator.isout, ode_integrator.reeval_fsal,
                               ode_integrator.u_modified, opts, destats, history,
                               ode_integrator)

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
