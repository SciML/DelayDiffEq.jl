"""
    savevalues!(integrator::DDEIntegrator, force_save=false)

Update solution of `integrator`, if necessary or forced by `force_save`.
"""
function savevalues!(integrator::DDEIntegrator, force_save=false)
    # update time of ODE integrator (can be slightly modified (< 10ϵ) because of time stops)
    if typeof(integrator.t) <: AbstractFloat # does not work for units!
        if integrator.integrator.t != integrator.t
            if abs(integrator.t - integrator.integrator.t) >= 10eps(integrator.t)
                error("unexpected time discrepancy detected")
            end

            integrator.integrator.t = integrator.t
            integrator.integrator.dt = integrator.integrator.t - integrator.integrator.tprev
        end
    end

    # update solution
    savevalues!(integrator.integrator, force_save, false) # reduce_size = false

    # update prev2_idx to indices of tprev and u(tprev) in solution
    # allows reset of ODE integrator (and hence history function) to the last
    # successful time step after failed steps
    integrator.prev2_idx = integrator.prev_idx

    # cache dt of interval [tprev, t] of ODE integrator since it can only be retrieved by
    # a possibly incorrect subtraction
    # NOTE: does not interfere with usual use of dtcache for non-adaptive methods since ODE
    # integrator is only used for inter- and extrapolation of future values and saving of
    # the solution but does not affect the size of time steps
    integrator.integrator.dtcache = integrator.integrator.dt

    # reduce ODE solution
    if !(typeof(integrator.saveat) <: Void)
        # obtain constant lags
        if typeof(integrator.sol.prob) <: ConstantLagDDEProblem
            #warn("ConstantLagDDEProblem is deprecated. Use DDEProblem instead.")
            constant_lags = integrator.sol.prob.lags
        else
            constant_lags = integrator.sol.prob.constant_lags
        end

        # delete part of ODE solution that is not required for DDE solution
        reduce_solution!(integrator,
                         # function values at later time points might be necessary for
                         # calculation of next step, thus keep those interpolation data
                         integrator.t - maximum(constant_lags))
    end

    # prevent reset of ODE integrator to cached values in the calculation of the next step
    # NOTE: does not interfere with usual use of accept_step since ODE integrator is only
    # used for inter- and extrapolation of future values and saving of the solution but does
    # not affect whether time steps are accepted
    integrator.integrator.accept_step = true
end

"""
    postamble!(integrator::DDEIntegrator)

Clean up solution of `integrator`.
"""
function postamble!(integrator::DDEIntegrator)
    # clean up solution of ODE integrator
    postamble!(integrator.integrator)

    # reduce solution if possible
    !(typeof(integrator.saveat) <: Void) && reduce_solution!(integrator,
                                                             integrator.sol.t[end])
end

"""
    perform_step!(integrator::DDEIntegrator)

Calculate next step of `integrator`.
"""
@muladd function perform_step!(integrator::DDEIntegrator)
    # reset ODE integrator to cached values if last step failed
    if !integrator.integrator.accept_step
        if typeof(integrator.u) <: AbstractArray
            recursivecopy!(integrator.integrator.u, integrator.sol.u[end])
        else
            integrator.integrator.u = integrator.sol.u[end]
        end
        integrator.integrator.t = integrator.sol.t[end]
        integrator.integrator.tprev = integrator.sol.t[integrator.prev2_idx]
        integrator.integrator.dt = integrator.integrator.dtcache

        # u(tprev) is not modified hence we do not have to copy it
        integrator.integrator.uprev = integrator.sol.u[integrator.prev2_idx]

        # do not have to reset interpolation data in initial time step since always a
        # constant extrapolation is used (and interpolation data of solution at initial
        # time point is not complete!)
        if length(integrator.sol.t) > 1
            recursivecopy!(integrator.integrator.k, integrator.sol.k[end])
        end
    end

    # reset boolean which indicates whether history function was evaluated at a time point
    # past the final point of the current solution
    # NOTE: does not interfere with usual use of isout since ODE integrator is only used for
    # inter- and extrapolation of future values and saving of the solution but does not
    # affect whether time steps are accepted
    integrator.integrator.isout = false

    # perform always at least one calculation
    perform_step!(integrator, integrator.cache)

    # if the history function was evaluated at time points past the final time point of the
    # solution, i.e. returned extrapolated values, continue with a fixed-point iteration
    if integrator.integrator.isout
        # update ODE integrator to next time interval together with correct interpolation
        advance_ode_integrator!(integrator)

        numiters = 1

        while true
            # calculate next step
            perform_step!(integrator, integrator.cache, true) # repeat_step=true

            # calculate residuals of fixed-point iteration
            # TODO: replace with updated calculate_residuals
            if typeof(integrator.u) <: AbstractArray
                @. integrator.resid = (integrator.u - integrator.integrator.u) /
                    (integrator.fixedpoint_abstol + max(abs(integrator.u),
                                                        abs(integrator.integrator.u)) *
                     integrator.fixedpoint_reltol)
            else
                integrator.resid = (integrator.u - integrator.integrator.u) /
                    (integrator.fixedpoint_abstol + max(abs(integrator.u),
                                                        abs(integrator.integrator.u)) *
                            integrator.fixedpoint_reltol)
            end

            # update error estimate of integrator with a combined error
            # estimate of both integrator and fixed-point iteration
            # this prevents acceptance of steps with poor performance in fixed-point
            # iteration
            integrator.EEst = max(integrator.EEst,
                                  integrator.fixedpoint_norm(integrator.resid))

            # complete interpolation data of DDE integrator for time interval [t, t+dt]
            # and copy it to ODE integrator
            # has to be done before updates to ODE integrator, otherwise history function
            # is incorrect
            if typeof(integrator.cache) <: OrdinaryDiffEq.CompositeCache
                OrdinaryDiffEq.ode_addsteps!(integrator.k, integrator.t, integrator.uprev,
                                             integrator.u, integrator.dt, integrator.f,
                                             integrator.cache.caches[integrator.cache.current],
                                             Val{false}, Val{true}, Val{true})
            else
                OrdinaryDiffEq.ode_addsteps!(integrator.k, integrator.t, integrator.uprev,
                                             integrator.u, integrator.dt, integrator.f,
                                             integrator.cache, Val{false}, Val{true},
                                             Val{true})
            end
            recursivecopy!(integrator.integrator.k, integrator.k)

            # update value u(t+dt)
            if typeof(integrator.u) <: AbstractArray
                recursivecopy!(integrator.integrator.u, integrator.u)
            else
                integrator.integrator.u = integrator.u
            end

            # stop fixed-point iteration when error estimate is small or maximal number of
            # steps is exceeded
            if integrator.EEst <= 1 || numiters > integrator.max_fixedpoint_iters
                break
            end

            numiters += 1
        end
    else
        # update ODE integrator to next time interval together with correct interpolation
        advance_ode_integrator!(integrator)
    end

    # is reset if time step is saved, prevents unnecessary copies and assignments
    integrator.integrator.accept_step = false
end

"""
    initialize!(integrator::DDEIntegrator)

Set initial values of `integrator`.
"""
function initialize!(integrator::DDEIntegrator)
    # initialize DDE integrator
    initialize!(integrator, integrator.cache)

    # copy interpolation data to ODE integrator
    integrator.integrator.kshortsize = integrator.kshortsize
    integrator.integrator.k = recursivecopy(integrator.k)

    # add interpolation steps to ODE integrator to ensure that interpolation data
    # is always maximal when calculating the next step
    # exact values do not matter since in the initial time step always a constant
    # extrapolation is used
    OrdinaryDiffEq.ode_addsteps!(integrator.integrator, integrator.f)
end

"""
    u_modified!(integrator::DDEIntegrator, bool::Bool)

Signal `integrator` whether state vector `u` was modified by a callback.

A modified `u` will lead to recalculations in order to prevent discontinuities.
"""
@inline function u_modified!(integrator::DDEIntegrator, bool::Bool)
    integrator.u_modified = bool
end

"""
    get_proposed_dt(integrator::DDEIntegrator)

Get the time step that `integrator` will take after the current step.
"""
@inline get_proposed_dt(integrator::DDEIntegrator) = integrator.dtpropose

"""
    set_proposed_dt!(integrator::DDEIntegrator, dt)

Set the time step that `integrator` will take after the current step to `dt`.
"""
@inline set_proposed_dt!(integrator::DDEIntegrator, dt) = (integrator.dtpropose = dt)

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

"""
    terminate!(integrator::DDEIntegrator)

Stop further calculations of `integrator`.
"""
function terminate!(integrator::DDEIntegrator)
    integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end
