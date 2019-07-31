"""
    advance_ode_integrator!(integrator::DDEIntegrator)

Advance the ODE integrator of `integrator` to the next time interval by updating its values
and interpolation data with the current values and a full set of interpolation data of
`integrator`.
"""
function advance_ode_integrator!(integrator::DDEIntegrator)
  @unpack f, u, t, p, k, dt, uprev, alg, cache = integrator
  ode_integrator = integrator.integrator

  # algorithm only works if current time of DDE integrator equals final time point
  # of solution
  t != ode_integrator.sol.t[end] && error("cannot advance ODE integrator")

  # complete interpolation data of DDE integrator for time interval [t, t+dt]
  # and copy it to ODE integrator
  # has to be done before updates to ODE integrator, otherwise history function
  # is incorrect
  if iscomposite(alg)
    DiffEqBase.addsteps!(k, t, uprev, u, dt, f, p, cache.caches[cache.current],
                         false, true, true)
  else
    DiffEqBase.addsteps!(k, t, uprev, u, dt, f, p, cache, false, true, true)
  end
  @inbounds for i in 1:length(k)
    copyat_or_push!(ode_integrator.k, i, k[i])
  end

  # move ODE integrator to interval [t, t+dt]
  ode_integrator.t = t + dt
  ode_integrator.tprev = t
  ode_integrator.dt = dt
  if iscomposite(alg)
    ode_integrator.cache.current = cache.current
  end
  if isinplace(integrator.sol.prob)
    recursivecopy!(ode_integrator.u, u)
  else
    ode_integrator.u = u
  end

  # u(t) is not modified hence we do not have to copy it
  ode_integrator.uprev = ode_integrator.sol.u[end]

  # update prev_idx to index of t and u(t) in solution
  integrator.prev_idx = length(ode_integrator.sol.t)

  nothing
end

"""
    move_back_ode_integrator!(integrator::DDEIntegrator)

Move the ODE integrator of `integrator` one integration step back by reverting its values
and interpolation data to the values saved in the dense history.
"""
function move_back_ode_integrator!(integrator::DDEIntegrator)
  ode_integrator = integrator.integrator
  @unpack sol = ode_integrator

  # set values of the ODE integrator back to the values in the solution
  if isinplace(sol.prob)
    recursivecopy!(ode_integrator.u, sol.u[end])
  else
    ode_integrator.u = sol.u[end]
  end
  ode_integrator.t = sol.t[end]
  ode_integrator.tprev = sol.t[integrator.prev2_idx]

  # u(tprev) is not modified hence we do not have to copy it
  ode_integrator.uprev = sol.u[integrator.prev2_idx]

  # revert to the previous time step
  ode_integrator.dt = ode_integrator.dtcache

  # we do not have to reset the interpolation data in the initial time step since always a
  # constant extrapolation is used (and interpolation data of solution at initial
  # time point is not complete!)
  if length(sol.t) > 1
    recursivecopy!(ode_integrator.k, sol.k[end])
  end

  nothing
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

# handle discontinuities at the current time point of the `integrator`
function OrdinaryDiffEq.handle_discontinuities!(integrator::DDEIntegrator)
    # remove all discontinuities at current time point and calculate minimal order
    # of these discontinuities
    d = pop!(integrator.opts.d_discontinuities)
    order = d.order
    while !isempty(integrator.opts.d_discontinuities) &&
        top(integrator.opts.d_discontinuities) == integrator.tdir * integrator.t

        d2 = pop!(integrator.opts.d_discontinuities)
        order = min(order, d2.order)
    end

    # remove all discontinuities close to the current time point as well and
    # calculate minimal order of these discontinuities
    # integrator.EEst has unitless type of integrator.t
    if typeof(integrator.EEst) <: AbstractFloat
        maxΔt = 10eps(integrator.t)

        while !isempty(integrator.opts.d_discontinuities) &&
            abs(top(integrator.opts.d_discontinuities).t - integrator.tdir * integrator.t) < maxΔt

            d2 = pop!(integrator.opts.d_discontinuities)
            order = min(order, d2.order)
        end

        # also remove all corresponding time stops
        while !isempty(integrator.opts.tstops) &&
            abs(top(integrator.opts.tstops) - integrator.tdir * integrator.t) < maxΔt

            pop!(integrator.opts.tstops)
        end
    end

    # add discontinuities of next order to integrator
    add_next_discontinuities!(integrator, order)

  nothing
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
  neutral = integrator.sol.prob.neutral
  next_order = neutral ? order : order + 1

  # only track discontinuities up to order of the applied method
  alg_maximum_order = OrdinaryDiffEq.alg_maximum_order(integrator.alg)
  next_order <= alg_maximum_order + 1 || return

  # discontinuities caused by constant lags
  if has_constant_lags(integrator)
    constant_lags = integrator.sol.prob.constant_lags
    maxlag = integrator.tdir * (integrator.sol.prob.tspan[end] - t)

    for lag in constant_lags
      if integrator.tdir * lag < maxlag
        # calculate discontinuity and add it to heap of discontinuities and time stops
        d = Discontinuity(integrator.tdir * (t + lag), next_order)
        push!(integrator.opts.d_discontinuities, d)
        push!(integrator.opts.tstops, d.t)
      end
    end
  end

  # track propagated discontinuities with callback
  push!(integrator.tracked_discontinuities, Discontinuity(integrator.tdir * t, order))

  nothing
end
