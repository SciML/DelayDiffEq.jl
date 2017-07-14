"""
    savevalues!(integrator::DDEIntegrator, force_save=false)

Update solution of `integrator`, if necessary or forced by `force_save`.
"""
function savevalues!(integrator::DDEIntegrator, force_save=false)
    # update ODE integrator
    integrator.integrator.u = integrator.u
    integrator.integrator.k = integrator.k
    integrator.integrator.t = integrator.t

    # update solution of DDE integrator
    invoke(savevalues!, Tuple{AbstractDDEIntegrator,Bool}, integrator, force_save)

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
    integrator.integrator.u = integrator.u
    integrator.integrator.k = integrator.k
    integrator.integrator.t = integrator.t

    # clean up solution of DDE integrator
    OrdinaryDiffEq.postamble!(integrator)

    # clean up solution of ODE integrator
    OrdinaryDiffEq.postamble!(integrator.integrator)
end

"""
    perform_step!(integrator::DDEIntegrator)

Calculate next step of `integrator`.
"""
function perform_step!(integrator::DDEIntegrator)
    # update previous time to extrapolate from current interval
    integrator.tprev = integrator.t

    # update ODE integrator
    integrator.integrator.uprev = integrator.uprev
    integrator.integrator.tprev = integrator.tprev
    integrator.integrator.fsalfirst = integrator.fsalfirst
    integrator.integrator.t = integrator.t
    integrator.integrator.dt = integrator.dt

    # if dt is greater than the minimal lag, then it's explicit so use fixed-point iteration
    if integrator.dt > minimum(integrator.prob.lags)

        # save these values to correct the extrapolation after the last iteration
        tprev_cache = integrator.tprev
        t_cache = integrator.t # same as tprev_cache?
        uprev_cache = integrator.uprev

        numiters = 1
        while true

            # save u to calculate residuals (u is overwritten in calculation of next step)
            if typeof(integrator.u) <: AbstractArray
                copy!(integrator.u_cache, integrator.u)
            else
                integrator.u_cache = integrator.u
            end

            # calculate next step
            perform_step!(integrator, integrator.cache)

            # calculate residuals
            if typeof(integrator.resid) <: AbstractArray
                @. integrator.resid = (integrator.u - integrator.u_cache) /
                    @muladd(integrator.fixedpoint_abstol + max(abs(integrator.u),
                                                               abs(integrator.u_cache)) *
                            integrator.fixedpoint_reltol)
            else
                integrator.resid = @. (integrator.u - integrator.u_cache) /
                    @muladd(integrator.fixedpoint_abstol + max(abs(integrator.u),
                                                               abs(integrator.u_cache)) *
                            integrator.fixedpoint_reltol)
            end

            # stop fixed-point iteration when residuals are small or maximal number of steps is exceeded
            fixedpointEEst = integrator.fixedpoint_norm(integrator.resid)
            if fixedpointEEst < 1 || numiters > integrator.max_fixedpoint_iters
                break
            end

            # special updates of ODE integrator after the first iteration step
            # to use interpolation of ODE integrator in the next iterations
            # when evaluating the history function
            if numiters == 1
                integrator.integrator.tprev = integrator.t
                integrator.integrator.t = integrator.t + integrator.dt
                integrator.integrator.uprev = integrator.u
            end

            # update ODE integrator
            integrator.integrator.u = integrator.u
            integrator.integrator.k = integrator.k

            numiters += 1
        end

        # reset values of DDE integrator after last iteration
        integrator.t = t_cache
        integrator.tprev = tprev_cache
        integrator.uprev = uprev_cache

        # update current time of ODE integrator
        integrator.integrator.t = t_cache
    else # no iterations
        perform_step!(integrator, integrator.cache)
    end

    #integrator.u = integrator.integrator.u
    #integrator.fsallast = integrator.integrator.fsallast
    #if integrator.opts.adaptive
    #  integrator.EEst = integrator.integrator.EEst
    #end
end

"""
    initialize!(integrator::DDEIntegrator)

Set initial values of `integrator`.
"""
function initialize!(integrator::DDEIntegrator)
    initialize!(integrator, integrator.cache, integrator.f)

    # set also initial values of ODE integrator
    initialize!(integrator.integrator, integrator.cache, integrator.f)
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
