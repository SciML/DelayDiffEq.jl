"""
    savevalues!(integrator::DDEIntegrator, force_save=false)

Update solution of `integrator`, if necessary or forced by `force_save`.
"""
function savevalues!(integrator::DDEIntegrator, force_save=false)
    # update ODE integrator to interval [tprev, t] with corresponding values
    # u(tprev) and u(t), and interpolation data k of this interval
    integrator.integrator.t = integrator.t
    integrator.integrator.tprev = integrator.tprev

    # copy u(tprev) since it is overwritten by integrator at the end of apply_step!
    if typeof(integrator.u) <: AbstractArray
        recursivecopy!(integrator.integrator.u, integrator.u)
        recursivecopy!(integrator.integrator.uprev, integrator.uprev)
    else
        integrator.integrator.u = integrator.u
        integrator.integrator.uprev = integrator.uprev
    end

    # copy interpolation data (fsalfirst overwritten at the end of apply_step!, which also
    # updates k[1] when using chaches for which k[1] points to fsalfirst)
    recursivecopy!(integrator.integrator.k, integrator.k)

    # add steps for interpolation to ODE integrator when needed
    OrdinaryDiffEq.ode_addsteps!(integrator.integrator, integrator.f)

    # update solution of ODE integrator
    savevalues!(integrator.integrator, force_save)
end

"""
    postamble!(integrator::DDEIntegrator)

Clean up solution of `integrator`.
"""
function postamble!(integrator::DDEIntegrator)
    # update ODE integrator
    integrator.integrator.t = integrator.t
    if typeof(integrator.u) <: AbstractArray
        recursivecopy!(integrator.integrator.u, integrator.u)
    else
        integrator.integrator.u = integrator.u
    end
    recursivecopy!(integrator.integrator.k, integrator.k)

    # add steps for interpolation to ODE integrator when needed
    OrdinaryDiffEq.ode_addsteps!(integrator.integrator, integrator.f)

    # clean up solution of ODE integrator
    OrdinaryDiffEq.postamble!(integrator.integrator)
end

"""
    perform_step!(integrator::DDEIntegrator)

Calculate next step of `integrator`.
"""
function perform_step!(integrator::DDEIntegrator)
    # cache error estimate of integrator and interpolation data of interval [tprev, t]
    # (maybe with already updated entry k[1] = fsalfirst == fsallast, if k[1] points to
    # fsalfirst) to be able to reset the corresponding variables in case calculation results
    # in numbers that are not finite
    recursivecopy!(integrator.k_cache, integrator.k)
    integrator.integrator.EEst = integrator.EEst

    # add steps to interpolation data ODE integrator if necessary
    OrdinaryDiffEq.ode_addsteps!(integrator.integrator, integrator.f)

    # perform always at least one calculation
    perform_step!(integrator, integrator.cache)

    # shrink interpolation data of ODE problem
    resize!(integrator.integrator.k, integrator.integrator.kshortsize)

    # if dt is greater than the minimal lag, then use a fixed-point iteration
    if integrator.dt > minimum(integrator.prob.lags) && isfinite(integrator.EEst)

        # update cached error estimate of integrator
        integrator.integrator.EEst = integrator.EEst

        # save value u(tprev) and interpolation data of interval [tprev, t] of ODE
        # integrator since they are overwritten by fixed-point iteration
        if typeof(integrator.uprev_cache) <: AbstractArray
            recursivecopy!(integrator.uprev_cache, integrator.integrator.uprev)
        else
            integrator.uprev_cache = integrator.integrator.uprev
        end
        recursivecopy!(integrator.k_integrator_cache, integrator.integrator.k)

        # move ODE integrator to interval [t, t+dt] to use interpolation of ODE integrator
        # in the next iterations when evaluating the history function
        integrator.integrator.t = integrator.t + integrator.dt
        integrator.integrator.tprev = integrator.t
        if typeof(integrator.integrator.uprev) <: AbstractArray
            recursivecopy!(integrator.integrator.uprev, integrator.uprev)
        else
            integrator.integrator.uprev = integrator.uprev
        end

        numiters=1

        while true
            # update value u(t+dt) and interpolation data of interval [t, t+dt] that are
            # used for the interpolation of the history function in the next iteration
            if typeof(integrator.u) <: AbstractArray
                recursivecopy!(integrator.integrator.u, integrator.u)
            else
                integrator.integrator.u = integrator.u
            end
            recursivecopy!(integrator.integrator.k, integrator.k)

            # add steps to interpolation data of ODE integrator if necessary
            OrdinaryDiffEq.ode_addsteps!(integrator.integrator, integrator.f)

            # calculate next step
            perform_step!(integrator, integrator.cache)

            # shrink interpolation data of ODE integrator
            resize!(integrator.integrator.k, integrator.integrator.kshortsize)

            # calculate residuals of fixed-point iteration
            if typeof(integrator.resid) <: AbstractArray
                @. integrator.resid = (integrator.u - integrator.integrator.u) /
                    @muladd(integrator.fixedpoint_abstol + max(abs(integrator.u),
                                                               abs(integrator.integrator.u)) *
                            integrator.fixedpoint_reltol)
            else
                integrator.resid = @. (integrator.u - integrator.integrator.u) /
                    @muladd(integrator.fixedpoint_abstol + max(abs(integrator.u),
                                                               abs(integrator.integrator.u)) *
                            integrator.fixedpoint_reltol)
            end
            fixedpointEEst = integrator.fixedpoint_norm(integrator.resid)

            # stop fixed-point iteration when error estimate of integrator or error estimate
            # of fixed-point iteration are not finite
            if !(isfinite(fixedpointEEst)) || !(isfinite(integrator.EEst))
                # assure that integrator is reset to cached values
                integrator.EEst = max(fixedpointEEst, integrator.EEst)
                break
            end

            # update cached value of error estimate of integrator with a combined error
            # estimate of both integrator and fixed-point iteration
            # this prevents acceptance of steps with poor performance in fixed-point
            # iteration
            integrator.integrator.EEst = max(fixedpointEEst, integrator.EEst)

            # stop fixed-point iteration when error estimate is small or maximal number of
            # steps is exceeded
            if integrator.integrator.EEst <= 1 || numiters > integrator.max_fixedpoint_iters
                # update error estimate with combined error estimate
                integrator.EEst = integrator.integrator.EEst
                break
            end

            numiters += 1
        end

        # reset ODE integrator to interval [tprev, t] with corresponding values
        # u(tprev) and u(t), and interpolation data k of this interval
        integrator.integrator.t = integrator.t
        integrator.integrator.tprev = integrator.tprev
        if typeof(integrator.u) <: AbstractArray
            recursivecopy!(integrator.integrator.u, integrator.uprev)
            recursivecopy!(integrator.integrator.uprev, integrator.uprev_cache)
        else
            integrator.integrator.u = integrator.uprev
            integrator.integrator.uprev = integrator.uprev_cache
        end
        recursivecopy!(integrator.integrator.k, integrator.k_integrator_cache)
    end

    # if error estimate of integrator is not a finite number reset it to last cached error
    # estimate or 2, and reset interpolation data of integrator to interpolation data of
    # interval [tprev, t] (maybe with updated entry k[1]) before calculation of current step
    # then current step will not be accepted, time step dt will be decreased,
    # and calculation of next step will be repeated starting with the same
    # initial interpolation data
    if !isfinite(integrator.EEst)
        integrator.EEst = max(2, integrator.integrator.EEst) # EEst must be > 1
        recursivecopy!(integrator.k, integrator.k_cache)
    end
end

"""
    initialize!(integrator::DDEIntegrator)

Set initial values of `integrator`.
"""
function initialize!(integrator::DDEIntegrator)
    initialize!(integrator, integrator.cache, integrator.f)

    # interpolation data of integrator and ODE integrator have to be cached
    # when next step is calculated
    integrator.k_cache = recursivecopy(integrator.k)
    integrator.k_integrator_cache = recursivecopy(integrator.k)

    # copy interpolation data to ODE integrator
    integrator.integrator.kshortsize = integrator.kshortsize
    integrator.integrator.k = recursivecopy(integrator.k)
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

Set the time step that `integrator` will take after the current step to `dt.
"""
@inline set_proposed_dt!(integrator::DDEIntegrator, dt) = (integrator.dtpropose = dt)

user_cache(integrator::DDEIntegrator) = user_cache(integrator)
u_cache(integrator::DDEIntegrator) = u_cache(integrator.cache)
du_cache(integrator::DDEIntegrator)= du_cache(integrator.cache)
full_cache(integrator::DDEIntegrator) = chain(u_cache(integrator), du_cache(integrator.cache))

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
