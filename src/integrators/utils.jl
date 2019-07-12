"""
    advance_ode_integrator!(integrator::DDEIntegrator)

Advance ODE integrator of `integrator` to next time interval, values and complete
interpolation data of `integrator`.
"""
function advance_ode_integrator!(integrator::DDEIntegrator)
    # algorithm only works if current time of DDE integrator equals final time point
    # of solution
    integrator.t != integrator.sol.t[end] && error("cannot advance ODE integrator")

    # complete interpolation data of DDE integrator for time interval [t, t+dt]
    # and copy it to ODE integrator
    # has to be done before updates to ODE integrator, otherwise history function
    # is incorrect
    if typeof(integrator.cache) <: OrdinaryDiffEq.CompositeCache
        addsteps!(integrator.k, integrator.t, integrator.uprev,
                                     integrator.u, integrator.dt, integrator.f,
                                     integrator.p,
                                     integrator.cache.caches[integrator.cache.current],
                                     false, true, true)
    else
        addsteps!(integrator.k, integrator.t, integrator.uprev,
                                     integrator.u, integrator.dt, integrator.f,
                                     integrator.p,
                                     integrator.cache, false, true,
                                     true)
    end
    recursivecopy!(integrator.integrator.k, integrator.k)

    # move ODE integrator to interval [t, t+dt]
    integrator.integrator.t = integrator.t + integrator.dt
    integrator.integrator.tprev = integrator.t
    integrator.integrator.dt = integrator.dt
    if typeof(integrator.alg) <: OrdinaryDiffEq.OrdinaryDiffEqCompositeAlgorithm
      integrator.integrator.cache.current = integrator.cache.current
    end
    if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.integrator.u, integrator.u)
    else
        integrator.integrator.u = integrator.u
    end

    # u(t) is not modified hence we do not have to copy it
    integrator.integrator.uprev = integrator.sol.u[end]

    # update prev_idx to index of t and u(t) in solution
    integrator.prev_idx = length(integrator.sol.t)
end

#=
Dealing with discontinuities

If we hit a discontinuity (this is checked in `apply_step!`), then we remove the
discontinuity, additional discontinuities at the current time point (if present), and
maybe also discontinuities and time stops coming shortly after the current time point
in `handle_discontinuities!`. The order of the discontinuity at the current time point is
defined as the lowest order of all these discontinuities.

If the problem is not neutral, we will only add additional discontinuities if
this order is less or equal to the order of the algorithm in
`add_next_discontinuities!`. If we add discontinuities, we add discontinuities
of the next order caused by constant lags (these we can calculate explicitly and
just add them to `d_discontinuities` and `tstops`) and we add the current
discontinuity to `tracked_discontinuities` which is the array of old
discontinuities that are checked by a `DiscontinuityCallback` (if existent).
=#

"""
    handle_discontinuities!(integrator::DDEIntegrator)

Handle discontinuities at the current time point of `integrator`.
"""
function handle_discontinuities!(integrator::DDEIntegrator)
    # remove all discontinuities at current time point and calculate minimal order
    # of these discontinuities
    d = pop!(integrator.opts.d_discontinuities)
    order = d.order
    while !isempty(integrator.opts.d_discontinuities) &&
        top(integrator.opts.d_discontinuities) == integrator.t

        d2 = pop!(integrator.opts.d_discontinuities)
        order = min(order, d2.order)
    end

    # remove all discontinuities close to the current time point as well and
    # calculate minimal order of these discontinuities
    # integrator.EEst has unitless type of integrator.t
    if typeof(integrator.EEst) <: AbstractFloat
        maxΔt = 10eps(integrator.t)

        while !isempty(integrator.opts.d_discontinuities) &&
            abs(top(integrator.opts.d_discontinuities).t - integrator.t) < maxΔt

            d2 = pop!(integrator.opts.d_discontinuities)
            order = min(order, d2.order)
        end

        # also remove all corresponding time stops
        while !isempty(integrator.opts.tstops) &&
            abs(top(integrator.opts.tstops) - integrator.t) < maxΔt

            pop!(integrator.opts.tstops)
        end
    end

    # add discontinuities of next order to integrator
    add_next_discontinuities!(integrator, order)
end

"""
    add_next_discontinuities!(integrator::DDEIntegrator, order[, t=integrator.t])

Add discontinuities of next order that are propagated from discontinuity of
order `order` at time `t` in `integrator`, but only if `order` is less or equal
than the order of the applied method or the problem is neutral.

Discontinuities caused by constant delays are immediately calculated, and
discontinuities caused by dependent delays are tracked by a callback.
"""
function add_next_discontinuities!(integrator, order, t=integrator.t)
    # obtain delays
    neutral = integrator.sol.prob.neutral
    constant_lags = integrator.sol.prob.constant_lags

    # only track discontinuities up to order of the applied method
    order > alg_maximum_order(integrator.alg) && !neutral && return nothing

    # discontinuities caused by constant lags
    if constant_lags !== nothing
        maxlag = abs(integrator.sol.prob.tspan[end] - t)
        for lag in constant_lags
            if abs(lag) < maxlag
                # calculate discontinuity and add it to heap of discontinuities and time stops
                d = Discontinuity(t + lag, order + 1)
                push!(integrator.opts.d_discontinuities, d)
                push!(integrator.opts.tstops, d.t)
            end
        end
    end

    # track propagated discontinuities with callback
    push!(integrator.tracked_discontinuities, Discontinuity(t, order))
end
