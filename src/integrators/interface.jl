# accept/reject computed integration step, propose next step, and apply callbacks
function OrdinaryDiffEq.loopfooter!(integrator::DDEIntegrator)
  # apply same logic as in OrdinaryDiffEq
  OrdinaryDiffEq._loopfooter!(integrator)

  if !integrator.accept_step
    # reset ODE integrator to the cached values if the last step failed
    move_back_ode_integrator!(integrator)

    # track propagated discontinuities for dependent delays
    if integrator.opts.adaptive && integrator.iter > 0 && has_dependent_lags(integrator)
      track_propagated_discontinuities!(integrator)
    end
  end

  nothing
end

# save current state of the integrator
function DiffEqBase.savevalues!(integrator::DDEIntegrator, force_save = false,
                                reduce_size = false)::Tuple{Bool,Bool}
  ode_integrator = integrator.integrator

  # update time of ODE integrator (can be slightly modified (< 10ϵ) because of time stops)
  # integrator.EEst has unitless type of integrator.t
  if typeof(integrator.EEst) <: AbstractFloat
    if ode_integrator.t != integrator.t
      abs(integrator.t - ode_integrator.t) < 10eps(integrator.t) ||
        error("unexpected time discrepancy detected")

      ode_integrator.t = integrator.t
      ode_integrator.dt = ode_integrator.t - ode_integrator.tprev
    end
  end

  # if forced, then the user or an event changed the state u directly.
  if force_save
    if isinplace(integrator.sol.prob)
      ode_integrator.u .= integrator.u
    else
      ode_integrator.u = integrator.u
    end
  end

  # update history
  saved, savedexactly = DiffEqBase.savevalues!(ode_integrator, force_save,
                                               false) # reduce_size = false

  # check that history was actually updated
  saved || error("dense history could not be updated")

  # update prev2_idx to indices of tprev and u(tprev) in solution
  # allows reset of ODE integrator (and hence history function) to the last
  # successful time step after failed steps
  integrator.prev2_idx = integrator.prev_idx

  # cache dt of interval [tprev, t] of ODE integrator since it can only be retrieved by
  # a possibly incorrect subtraction
  # NOTE: does not interfere with usual use of dtcache for non-adaptive methods since ODE
  # integrator is only used for inter- and extrapolation of future values and saving of
  # the solution but does not affect the size of time steps
  ode_integrator.dtcache = ode_integrator.dt

  # update solution
  OrdinaryDiffEq._savevalues!(integrator, force_save, reduce_size)
end

# clean up the solution of the integrator
function DiffEqBase.postamble!(integrator::DDEIntegrator)
  # clean up solution of the ODE integrator
  DiffEqBase.postamble!(integrator.integrator)

  # clean solution of the DDE integrator
  OrdinaryDiffEq._postamble!(integrator)
end

# perform next integration step
function OrdinaryDiffEq.perform_step!(integrator::DDEIntegrator)
  @unpack f, t, p, k, uprev, dt, resid, alg, cache, max_fixedpoint_iters = integrator
  @unpack fixedpoint_abstol, fixedpoint_reltol, fixedpoint_norm = integrator
  ode_integrator = integrator.integrator
  internalnorm = integrator.opts.internalnorm
  prob = integrator.sol.prob

  # reset boolean which indicates whether history function was evaluated at a time point
  # past the final point of the current solution
  # NOTE: does not interfere with usual use of isout since ODE integrator is only used for
  # inter- and extrapolation of future values and saving of the solution but does not
  # affect whether time steps are accepted
  ode_integrator.isout = false

  # perform always at least one calculation of the stages
  OrdinaryDiffEq.perform_step!(integrator, cache)

  # if the history function was evaluated at time points past the final time point of the
  # solution, i.e. returned extrapolated values, continue with a fixed-point iteration
  if ode_integrator.isout
    # update ODE integrator to next time interval together with correct interpolation
    advance_ode_integrator!(integrator)

    numiters = 1

    while true
      # calculate next step
      OrdinaryDiffEq.perform_step!(integrator, cache, true) # repeat_step=true

      # calculate residuals of fixed-point iteration
      if isinplace(prob)
        OrdinaryDiffEq.calculate_residuals!(resid, ode_integrator.u, integrator.u,
                                            fixedpoint_abstol, fixedpoint_reltol,
                                            internalnorm, t)
      else
        resid = OrdinaryDiffEq.calculate_residuals(ode_integrator.u, integrator.u,
                                                   fixedpoint_abstol, fixedpoint_reltol,
                                                   internalnorm, t)
      end

      # update error estimate of integrator with a combined error
      # estimate of both integrator and fixed-point iteration
      # this prevents acceptance of steps with poor performance in fixed-point
      # iteration
      integrator.EEst = max(integrator.EEst, fixedpoint_norm(resid, t))

      # complete interpolation data of DDE integrator for time interval [t, t+dt]
      # and copy it to ODE integrator
      # has to be done before updates to ODE integrator, otherwise history function
      # is incorrect
      if iscomposite(alg)
        DiffEqBase.addsteps!(k, t, uprev, integrator.u, dt, f, p,
                             cache.caches[cache.current], false, true, true)
      else
        DiffEqBase.addsteps!(k, t, uprev, integrator.u, dt, f, p, cache, false, true, true)
      end
      recursivecopy!(ode_integrator.k, k)

      # update value u(t+dt)
      if isinplace(prob)
        recursivecopy!(ode_integrator.u, integrator.u)
      else
        ode_integrator.u = integrator.u
      end

      # stop fixed-point iteration when error estimate is small or maximal number of
      # steps is exceeded
      if integrator.EEst <= 1 || numiters > max_fixedpoint_iters
        break
      end

      numiters += 1
    end
  else
    # update ODE integrator to next time interval together with correct interpolation
    advance_ode_integrator!(integrator)
  end

  nothing
end

# initialize the integrator
function OrdinaryDiffEq.initialize!(integrator::DDEIntegrator)
  ode_integrator = integrator.integrator

  # initialize the cache
  OrdinaryDiffEq.initialize!(integrator, integrator.cache)

  # copy interpolation data to the ODE integrator
  ode_integrator.kshortsize = integrator.kshortsize
  ode_integrator.k = recursivecopy(integrator.k)

  # add interpolation steps to ODE integrator to ensure that interpolation data
  # is always maximal when calculating the next step
  # exact values do not matter since in the initial time step always a constant
  # extrapolation is used
  DiffEqBase.addsteps!(ode_integrator, integrator.f)

  nothing
end

# signal the integrator that state u was modified
function DiffEqBase.u_modified!(integrator::DDEIntegrator, bool::Bool)
  integrator.u_modified = bool
end

# return next integration time step
DiffEqBase.get_proposed_dt(integrator::DDEIntegrator) = integrator.dtpropose

# set next integration time step
function DiffEqBase.set_proposed_dt!(integrator::DDEIntegrator, dt)
  integrator.dtpropose = dt
  nothing
end

DiffEqBase.get_tmp_cache(integrator::DDEIntegrator) =
  DiffEqBase.get_tmp_cache(integrator, integrator.alg, integrator.cache)
user_cache(integrator::DDEIntegrator) = user_cache(integrator.cache)
u_cache(integrator::DDEIntegrator) = u_cache(integrator.cache)
du_cache(integrator::DDEIntegrator)= du_cache(integrator.cache)
full_cache(integrator::DDEIntegrator) = chain(user_cache(integrator),u_cache(integrator),du_cache(integrator))

resize!(integrator::DDEIntegrator, i::Int) = resize!(integrator, integrator.cache, i)
function resize!(integrator::DDEIntegrator, cache, i)
    for c in full_cache(integrator)
        resize!(c, i)
    end
end

function resize!(integrator::DDEIntegrator, cache::Union{Rosenbrock23Cache,
                                                         Rosenbrock32Cache}, i)
    for c in full_cache(integrator)
        resize!(c, i)
    end
    for c in vecu_cache(integrator.cache)
        resize!(c, i)
    end
    Jvec = vec(cache.J)
    cache.J = reshape(resize!(Jvec, i*i), i, i)
    Wvec = vec(cache.W)
    cache.W = reshape(resize!(Wvec, i*i), i, i)
end

function resize!(integrator::DDEIntegrator, cache::Union{ImplicitEulerCache,TrapezoidCache},
                 i)
    for c in full_cache(integrator)
        resize!(c, i)
    end
    for c in vecu_cache(integrator.cache)
        resize!(c, i)
    end
    for c in dual_cache(integrator.cache)
        resize!(c.du, i)
        resize!(c.dual_du, i)
    end
    if alg_autodiff(integrator.alg)
        cache.adf = autodiff_setup(cache.rhs, cache.uhold, integrator.alg)
    end
end

function deleteat!(integrator::DDEIntegrator, i::Int)
    for c in full_cache(integrator)
        deleteat!(c, i)
    end
end

# terminate integration
function DiffEqBase.terminate!(integrator::DDEIntegrator, retcode = :Terminated)
  integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, retcode)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
  nothing
end

# integrator can be reinitialized
DiffEqBase.has_reinit(integrator::DDEIntegrator) = true

"""
    reinit!(integrator::DDEIntegrator[, u0 = integrator.sol.prob.u0;
            t0 = integrator.sol.prob.tspan[1],
            tf = integrator.sol.prob.tspan[2],
            erase_sol = true,
            kwargs...])

Reinitialize `integrator` with (optionally) different initial state `u0`, different
integration interval from `t0` to `tf`, and erased solution if `erase_sol = true`.
"""
function DiffEqBase.reinit!(integrator::DDEIntegrator, u0 = integrator.sol.prob.u0;
                            t0 = integrator.sol.prob.tspan[1],
                            tf = integrator.sol.prob.tspan[end],
                            erase_sol = true,
                            tstops = integrator.opts.tstops_cache,
                            saveat = integrator.opts.saveat_cache,
                            d_discontinuities = integrator.opts.d_discontinuities_cache,
                            order_discontinuity_t0 = t0 == integrator.sol.prob.tspan[1] && u0 == integrator.sol.prob.u0 ? integrator.order_discontinuity_t0 : 0,
                            reset_dt = iszero(integrator.dtcache) && integrator.opts.adaptive,
                            reinit_callbacks = true, initialize_save = true,
                            reinit_cache = true)
  # reinit history
  reinit!(integrator.integrator, u0;
          t0 = t0, tf = tf, erase_sol = true, reset_dt = false, reinit_callbacks = false,
          reinit_cache = false)
  integrator.integrator.dt = zero(integrator.dt)
  integrator.integrator.dtcache = zero(integrator.dt)

  # reinit initial values of the integrator
  if isinplace(integrator.sol.prob)
    recursivecopy!(integrator.u, u0)
    recursivecopy!(integrator.uprev, integrator.u)
  else
    integrator.u = u0
    integrator.uprev = integrator.u
  end

  if OrdinaryDiffEq.alg_extrapolates(integrator.alg)
    if isinplace(integrator.sol.prob)
      recursivecopy!(integrator.uprev2, integrator.uprev)
    else
      integrator.uprev2 = integrator.uprev
    end
  end
  integrator.t = t0
  integrator.tprev = t0

  # reinit time stops, time points at which solution is saved, and discontinuities
  maximum_order = OrdinaryDiffEq.alg_maximum_order(integrator.alg)
  tstops_internal, saveat_internal, d_discontinuities_internal =
    OrdinaryDiffEq.tstop_saveat_disc_handling(tstops, saveat, d_discontinuities,
                                              integrator.tdir, (t0, tf),
                                              order_discontinuity_t0, maximum_order,
                                              integrator.sol.prob.constant_lags,
                                              typeof(integrator.t))

  integrator.opts.tstops = tstops_internal
  integrator.opts.saveat = saveat_internal
  integrator.opts.d_discontinuities = d_discontinuities_internal

  # update order of initial discontinuity
  integrator.order_discontinuity_t0 = order_discontinuity_t0

  # erase solution
  if erase_sol
    # resize vectors in solution
    resize_start = integrator.opts.save_start ? 1 : 0
    resize!(integrator.sol.u, resize_start)
    resize!(integrator.sol.t, resize_start)
    resize!(integrator.sol.k, resize_start)
    iscomposite(integrator.alg) && resize!(integrator.sol.alg_choice, resize_start)
    integrator.sol.u_analytic !== nothing && resize!(integrator.sol.u_analytic, 0)

    # save initial values
    if integrator.opts.save_start
      copyat_or_push!(integrator.sol.t, 1, integrator.t)
      if integrator.opts.save_idxs === nothing
        copyat_or_push!(integrator.sol.u, 1, integrator.u)
      else
        u_initial = integrator.u[integrator.opts.save_idxs]
        copyat_or_push!(integrator.sol.u, 1, u_initial, Val{false})
      end
    end

    # reset iteration counter
    integrator.saveiter = resize_start

    # erase array of tracked discontinuities
    if order_discontinuity_t0 ≤ OrdinaryDiffEq.alg_maximum_order(integrator.alg)
      resize!(integrator.tracked_discontinuities, 1)
      integrator.tracked_discontinuities[1] = Discontinuity(integrator.t, order_discontinuity_t0)
    else
      resize!(integrator.tracked_discontinuities, 0)
    end

    # reset history counters
    integrator.prev_idx = 1
    integrator.prev2_idx = 1
  end

  # reset integration counters
  integrator.iter = 0
  integrator.success_iter = 0

  # full re-initialize the PI in timestepping
  integrator.qold = integrator.opts.qoldinit
  integrator.q11 = one(integrator.t)
  integrator.erracc = one(integrator.erracc)
  integrator.dtacc = one(integrator.dtacc)
  integrator.u_modified = false

  if reset_dt
    DiffEqBase.auto_dt_reset!(integrator)
  end

  if reinit_callbacks
    OrdinaryDiffEq.initialize_callbacks!(integrator, initialize_save)
  end

  if reinit_cache
    OrdinaryDiffEq.initialize!(integrator)
  end

  nothing
end













function DiffEqBase.auto_dt_reset!(integrator::DDEIntegrator)
  @unpack f, u, t, tdir, opts, sol = integrator
  @unpack prob = sol
  @unpack abstol, reltol, internalnorm = opts

  # determine maximal time step
  if has_constant_lags(prob)
    dtmax = tdir * min(abs(opts.dtmax), minimum(abs, prob.constant_lags))
  else
    dtmax = opts.dtmax
  end

  # determine initial time step
  ode_prob = ODEProblem(f, prob.u0, prob.tspan, prob.p)
  integrator.dt = OrdinaryDiffEq.ode_determine_initdt(
    u, t, tdir, dtmax, opts.abstol, opts.reltol, opts.internalnorm, ode_prob, integrator)

  nothing
end

function DiffEqBase.add_tstop!(integrator::DDEIntegrator,t)
  integrator.tdir * (t - integrator.t) < 0 && error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.tstops, t)
end

function DiffEqBase.add_saveat!(integrator::DDEIntegrator,t)
  integrator.tdir * (t - integrator.t) < 0 && error("Tried to add a saveat that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.saveat, t)
end

@inline function DiffEqBase.get_du(integrator::DDEIntegrator)
  integrator.fsallast
end

@inline function DiffEqBase.get_du!(out,integrator::DDEIntegrator)
  out .= integrator.fsallast
end

DiffEqBase.addsteps!(integrator::DDEIntegrator,args...) = OrdinaryDiffEq._ode_addsteps!(integrator,args...)
DiffEqBase.change_t_via_interpolation!(integrator::DDEIntegrator,
                                        t,modify_save_endpoint::Type{Val{T}}=Val{false}) where T =
                                          OrdinaryDiffEq._change_t_via_interpolation!(integrator,t,modify_save_endpoint)
