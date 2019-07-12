"""
    discontinuity_function(integrator::DDEIntegrator, lag, T, t)

Evaluate function ``f(x) = T + lag(u(x), p, x) - x`` at time point `t`, where `T` is time
point of a previous discontinuity and `lag` is a dependent delay.
"""
function discontinuity_function(integrator::DDEIntegrator, lag, T, t)
  tmp = get_tmp_cache(integrator)
  cache = tmp === nothing ? nothing : first(tmp)

  # estimate state at the given time
  if cache === nothing
    ut = integrator(t, Val{0})
  else
    integrator(cache, t, Val{0})
    ut = cache
  end

  T + lag(ut, integrator.p, t) - t
end

"""
    discontinuity_interval(integrator::DDEIntegrator, lag, T, Θs)

Return an estimated subinterval of the current integration step of the `integrator` that
contains a propagated discontinuity induced by the dependent delay `lag` and the
discontinuity at time point `T`, or `nothing`.

The interval is estimated by checking the signs of `T + lag(u(t), p, t) - t` for time points
`integrator.t .+ θs` in the interval `[integrator.t, integrator.t + integrator.dt]`.
"""
function discontinuity_interval(integrator::DDEIntegrator, lag, T, Θs)
  # use start and end point of last time interval to check for discontinuities
  previous_condition = discontinuity_function(integrator, lag, T, integrator.t)
  if isapprox(previous_condition, 0;
              rtol = integrator.discontinuity_reltol,
              atol = integrator.discontinuity_abstol)
    prev_sign = 0
  else
    prev_sign = cmp(previous_condition, zero(previous_condition))
  end
  new_condition = discontinuity_function(integrator, lag, T, integrator.t + integrator.dt)
  new_sign = cmp(new_condition, zero(new_condition))

  # if signs are different we already know that a discontinuity exists
  if prev_sign * new_sign < 0
    return (zero(eltype(Θs)), one(eltype(Θs)))
  end

  # recheck interpolation intervals if no discontinuity found yet
  prev_Θ = zero(eltype(Θs))

  for i in 2:length(Θs)
    # evaluate sign at next interpolation point
    new_Θ = Θs[i]
    new_t = integrator.t + new_Θ * integrator.dt
    new_condition = discontinuity_function(integrator, lag, T, new_t)
    new_sign = cmp(new_condition, zero(new_condition))

    # return estimated interval if we find a root or observe a switch of signs
    if new_sign == 0
      return (new_Θ, new_Θ)
    elseif prev_sign * new_sign < 0
      return (prev_Θ, new_Θ)
    else
      # otherwise update lower estimate of subinterval
      prev_sign = new_sign
      prev_Θ = new_Θ
    end
  end

  nothing
end

"""
    discontinuity_time(integrator::DDEIntegrator, lag, T, interval)

Estimate time point of the propagated discontinuity induced by the dependent delay
`lag` and the discontinuity at time point `T` inside the `interval` of the current
integration step of the `integrator`.
"""
function discontinuity_time(integrator::DDEIntegrator, lag, T, (bottom_Θ, top_Θ))
  if bottom_Θ == top_Θ
    # in that case we have already found the time point of a discontinuity
    Θ = top_Θ
  else
    # define function for root finding
    zero_func(Θ) = discontinuity_function(integrator, lag, T,
                                          integrator.t + Θ * integrator.dt)

    Θ = prevfloat(find_zero(zero_func, (bottom_Θ, top_Θ), Roots.AlefeldPotraShi();
                            atol = integrator.discontinuity_abstol / 100))
  end

  # Θ = prevfloat(...)
  # prevfloat guerentees that the new time is either 1 floating point
  # numbers just before the event or directly at zero, but not after.
  # If there's a barrier which is never supposed to be crossed,
  # then this will ensure that
  # The item never leaves the domain. Otherwise Roots.jl can return
  # a float which is slightly after, making it out of the domain, causing
  # havoc.

  integrator.t + Θ * integrator.dt
end


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
