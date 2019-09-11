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
function DiffEqBase.savevalues!(integrator::HistoryODEIntegrator, force_save = false,
                                reduce_size = false)::Tuple{Bool,Bool}
  integrator.saveiter += 1
  # TODO: save history only for a subset of components
  copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
  copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)

  integrator.saveiter_dense += 1
  copyat_or_push!(integrator.sol.k, integrator.saveiter_dense, integrator.k)


  if iscomposite(integrator.alg)
    copyat_or_push!(integrator.sol.alg_choice, integrator.saveiter, integrator.cache.current)
  end

  true, true
end

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
function DiffEqBase.postamble!(integrator::HistoryODEIntegrator)
  if integrator.saveiter == 0 || integrator.sol.t[integrator.saveiter] != integrator.t
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.t, integrator.saveiter, integrator.t)
    copyat_or_push!(integrator.sol.u, integrator.saveiter, integrator.u)

    integrator.saveiter_dense +=1
    copyat_or_push!(integrator.sol.k, integrator.saveiter_dense, integrator.k)

    if iscomposite(integrator.alg)
      copyat_or_push!(integrator.sol.alg_choice, integrator.saveiter, integrator.cache.current)
    end
  end

  resize!(integrator.sol.t, integrator.saveiter)
  resize!(integrator.sol.u, integrator.saveiter)
  resize!(integrator.sol.k, integrator.saveiter_dense)

  nothing
end

function DiffEqBase.postamble!(integrator::DDEIntegrator)
  # clean up solution of the ODE integrator
  DiffEqBase.postamble!(integrator.integrator)

  # clean solution of the DDE integrator
  OrdinaryDiffEq._postamble!(integrator)
end

# perform next integration step
function OrdinaryDiffEq.perform_step!(integrator::DDEIntegrator)
  @unpack cache, history = integrator

  # reset boolean which indicates if the history function was evaluated at a time point
  # past the final point of the current solution
  history.isout = false

  # perform always at least one calculation of the stages
  OrdinaryDiffEq.perform_step!(integrator, cache)

  # if the history function was evaluated at time points past the final time point of the
  # solution, i.e. returned extrapolated values, continue with a fixed-point iteration
  if history.isout
    # perform fixed-point iteration
    fpsolve!(integrator.fpsolver, integrator)
  end

  # update ODE integrator to next time interval together with correct interpolation
  advance_or_update_ode_integrator!(integrator)

  nothing
end

# initialize the integrator
function OrdinaryDiffEq.initialize!(integrator::DDEIntegrator)
  ode_integrator = integrator.integrator

  # initialize the cache
  OrdinaryDiffEq.initialize!(integrator, integrator.cache)

  # copy interpolation data to the ODE integrator
  @inbounds for i in 1:length(integrator.k)
    copyat_or_push!(ode_integrator.k, i, integrator.k[i])
  end

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

# obtain caches
DiffEqBase.get_tmp_cache(integrator::DDEIntegrator) =
  get_tmp_cache(integrator, integrator.alg, integrator.cache)
DiffEqBase.user_cache(integrator::DDEIntegrator) = user_cache(integrator.cache)
DiffEqBase.u_cache(integrator::DDEIntegrator) = u_cache(integrator.cache)
DiffEqBase.du_cache(integrator::DDEIntegrator) = du_cache(integrator.cache)
DiffEqBase.full_cache(integrator::DDEIntegrator) = full_cache(integrator.cache)

# change number of components
Base.resize!(integrator::DDEIntegrator, i::Int) = resize!(integrator, integrator.cache, i)
function Base.resize!(integrator::DDEIntegrator, cache, i)
  # resize ODE integrator (do only have to care about u and k)
  ode_integrator = integrator.integrator
  resize!(ode_integrator.u, i)
  for k in ode_integrator.k
    resize!(k, i)
  end

  # resize DDE integrator
  for c in full_cache(cache)
    resize!(c, i)
  end
  OrdinaryDiffEq.resize_nlsolver!(integrator, i)
  OrdinaryDiffEq.resize_J_W!(cache, integrator, i)
  resize_non_user_cache!(integrator, cache, i)
  fpsolve_resize!(integrator, i)
  nothing
end

DiffEqBase.resize_non_user_cache!(integrator::DDEIntegrator, cache, i) = nothing

function DiffEqBase.resize_non_user_cache!(integrator::DDEIntegrator,
                                           cache::RosenbrockMutableCache, i)
  cache.J = similar(cache.J, i, i)
  cache.W = similar(cache.W, i, i)
  cache.jac_config = DiffEqBase.resize_jac_config!(cache.jac_config, i)
  cache.grad_config = OrdinaryDiffEq.resize_grad_config!(cache.grad_config, i)
  nothing
end

function DiffEqBase.resize_non_user_cache!(integrator::DDEIntegrator,
                                           cache::Union{GenericImplicitEulerCache,GenericTrapezoidCache},
                                           i)
  cache.nl_rhs = integrator.alg.nlsolve(Val{:init}, cache.rhs, cache.u)
  nothing
end

# delete component(s)
function Base.deleteat!(integrator::DDEIntegrator, idxs)
  # delete components of ODE integrator (do only have to care about u and k)
  ode_integrator = integrator.integrator
  deleteat!(ode_integrator.u, idxs)
  for k in ode_integrator.k
    deleteat!(k, idxs)
  end

  # delete components of DDE integrator
  for c in full_cache(integrator)
    deleteat!(c, idxs)
  end
  deleteat_non_user_cache!(integrator, integrator.cache, i)
end

function DiffEqBase.deleteat_non_user_cache!(integrator::DDEIntegrator, cache, idxs)
  i = length(integrator.u)
  resize_non_user_cache!(integrator, cache, i)
end

# add component(s)
function DiffEqBase.addat!(integrator::DDEIntegrator, idxs)
  # add components to ODE integrator (do only have to care about u and k)
  ode_integrator = integrator.integrator
  addat!(ode_integrator.u, idxs)
  for k in ode_integrator.k
    addat!(k, idxs)
  end

  # add components to DDE integrator
  for c in full_cache(integrator)
    addat!(c, idxs)
  end
  addat_non_user_cache!(integrator, integrator.cache, idxs)
end

function DiffEqBase.addat_non_user_cache!(integrator::DDEIntegrator, cache, idxs)
  i = length(integrator.u)
  resize_non_user_cache!(integrator, cache, i)
end

# error check
DiffEqBase.last_step_failed(integrator::DDEIntegrator) =
  integrator.last_stepfail && !integrator.opts.adaptive

# terminate integration
function DiffEqBase.terminate!(integrator::DDEIntegrator, retcode = :Terminated)
  integrator.sol = DiffEqBase.solution_new_retcode(integrator.sol, retcode)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
  nothing
end

# integrator can be reinitialized
DiffEqBase.has_reinit(::HistoryODEIntegrator) = true
DiffEqBase.has_reinit(integrator::DDEIntegrator) = true

function DiffEqBase.reinit!(integrator::HistoryODEIntegrator, u0 = integrator.sol.prob.u0;
                            t0 = integrator.sol.prob.tspan[1],
                            tf = integrator.sol.prob.tspan[end],
                            erase_sol = true)
  # reinit initial values of the integrator
  if isinplace(integrator.sol.prob)
    recursivecopy!(integrator.u, u0)
    recursivecopy!(integrator.uprev, integrator.u)
  else
    integrator.u = u0
    integrator.uprev = integrator.u
  end
  integrator.t = t0
  integrator.tprev = t0
  integrator.dt = zero(integrator.dt)
  integrator.dtcache = zero(integrator.dtcache)

  # erase solution
  if erase_sol
    # resize vectors in solution
    resize!(integrator.sol.u, 1)
    resize!(integrator.sol.t, 1)
    resize!(integrator.sol.k, 1)
    iscomposite(integrator.alg) && resize!(integrator.sol.alg_choice, 1)

    # save initial values
    copyat_or_push!(integrator.sol.t, 1, integrator.t)
    copyat_or_push!(integrator.sol.u, 1, integrator.u)

    # reset iteration counter
    integrator.saveiter = 1
  end

  nothing
end

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
  reinit!(integrator.integrator, u0; t0 = t0, tf = tf, erase_sol = true)

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
  tType = typeof(integrator.t)
  tstops_internal, saveat_internal, d_discontinuities_internal =
    OrdinaryDiffEq.tstop_saveat_disc_handling(tstops, saveat, d_discontinuities,
                                              (tType(t0), tType(tf)),
                                              order_discontinuity_t0, maximum_order,
                                              integrator.sol.prob.constant_lags,
                                              integrator.sol.prob.neutral)

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
      integrator.tracked_discontinuities[1] = Discontinuity(integrator.tdir * integrator.t, order_discontinuity_t0)
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
  @unpack f, u, t, tdir, opts, sol, destats = integrator
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

  # update statistics
  destats.nf += 2

  nothing
end

function DiffEqBase.add_tstop!(integrator::DDEIntegrator,t)
  integrator.tdir * (t - integrator.t) < 0 && error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.tstops, integrator.tdir * t)
end

function DiffEqBase.add_saveat!(integrator::DDEIntegrator,t)
  integrator.tdir * (t - integrator.t) < 0 && error("Tried to add a saveat that is behind the current time. This is strictly forbidden")
  push!(integrator.opts.saveat, integrator.tdir * t)
end

@inline function DiffEqBase.get_du(integrator::DDEIntegrator)
  integrator.fsallast
end

@inline function DiffEqBase.get_du!(out,integrator::DDEIntegrator)
  out .= integrator.fsallast
end

DiffEqBase.has_destats(::DDEIntegrator) = true

DiffEqBase.addsteps!(integrator::DDEIntegrator,args...) = OrdinaryDiffEq._ode_addsteps!(integrator,args...)
DiffEqBase.change_t_via_interpolation!(integrator::DDEIntegrator,
                                        t,modify_save_endpoint::Type{Val{T}}=Val{false}) where T =
                                          OrdinaryDiffEq._change_t_via_interpolation!(integrator,t,modify_save_endpoint)

# update integrator when u is modified by callbacks
function OrdinaryDiffEq.handle_callback_modifiers!(integrator::DDEIntegrator)
    integrator.reeval_fsal = true # recalculate fsalfirst after applying step

    # update heap of discontinuities
    # discontinuity is assumed to be of order 0, i.e. solution x is discontinuous
    push!(integrator.opts.d_discontinuities, Discontinuity(integrator.tdir * integrator.t, 0))
end

# recalculate interpolation data and update the ODE integrator
function DiffEqBase.reeval_internals_due_to_modification!(
  integrator::DDEIntegrator, x::Type{Val{not_initialization}} = Val{true}) where not_initialization

  ode_integrator = integrator.integrator

  if not_initialization
    # update interpolation data of the integrator using the old dense history
    # of the ODE integrator
    resize!(integrator.k, integrator.kshortsize)
    DiffEqBase.addsteps!(integrator, integrator.f, true, true, true)

    # copy interpolation data to the ODE integrator
    recursivecopy!(ode_integrator.k, integrator.k)
  end

  # adjust current interval of the ODE integrator
  ode_integrator.t = integrator.t
  ode_integrator.dt = integrator.dt
  if isinplace(integrator.sol.prob)
    recursivecopy!(ode_integrator.u, integrator.u)
  else
    ode_integrator.u = integrator.u
  end

  integrator.u_modified = false
end

# perform one step
function DiffEqBase.step!(integrator::DDEIntegrator)
  @inbounds begin
    if integrator.opts.advance_to_tstop
      while integrator.tdir * integrator.t < top(integrator.opts.tstops)
        OrdinaryDiffEq.loopheader!(integrator)
        DiffEqBase.check_error!(integrator) === :Success || return
        OrdinaryDiffEq.perform_step!(integrator)
        OrdinaryDiffEq.loopfooter!(integrator)
      end
    else
      OrdinaryDiffEq.loopheader!(integrator)
      DiffEqBase.check_error!(integrator) === :Success || return
      OrdinaryDiffEq.perform_step!(integrator)
      OrdinaryDiffEq.loopfooter!(integrator)

      while !integrator.accept_step
        OrdinaryDiffEq.loopheader!(integrator)
        OrdinaryDiffEq.perform_step!(integrator)
        OrdinaryDiffEq.loopfooter!(integrator)
      end
    end
    OrdinaryDiffEq.handle_tstop!(integrator)
  end
end
