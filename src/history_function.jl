"""
    HistoryFunction(h, integrator)

Wrap history function `h` and integrator `integrator` to create a common
interface for retrieving values at any time point with varying accuracy.

Before the initial time point of the solution of the `integrator` values are calculated by
history function `h`, for time points up to the final time point of the solution
interpolated values of the solution are returned, and after the final time point an inter-
or extrapolation of the current state of the `integrator` is retrieved.
"""
mutable struct HistoryFunction{H,I<:ODEIntegrator} <: Function
  h::H
  integrator::I
  isout::Bool
end

HistoryFunction(h, integrator) = HistoryFunction(h, integrator, false)

function (f::HistoryFunction)(p, t, ::Type{Val{deriv}}=Val{0}; idxs=nothing) where deriv
  @unpack integrator = f
  @unpack sol = integrator

  @inbounds if integrator.tdir * t < integrator.tdir * sol.t[1]
    if deriv == 0 && idxs === nothing
      return f.h(p, t)
    elseif idxs === nothing
      return f.h(p, t, Val{deriv})
    elseif deriv == 0
      return f.h(p, t; idxs = idxs)
    else
      return f.h(p, t, Val{deriv}; idxs = idxs)
    end
  elseif integrator.tdir * t <= integrator.tdir * sol.t[end] # Put equals back
    return sol.interp(t, idxs, Val{deriv}, f.integrator.p)
  end

  # history function is evaluated at time point past the final time point of
  # the current solution
  f.isout = true

  # handle extrapolations at initial time point
  if integrator.t == sol.prob.tspan[1]
    return constant_extrapolant(t, integrator, idxs, Val{deriv})
  else
    return OrdinaryDiffEq.current_interpolant(t, integrator, idxs, Val{deriv})
  end
end

function (f::HistoryFunction)(val, p, t, ::Type{Val{deriv}}=Val{0}; idxs=nothing) where deriv
  @unpack integrator = f
  @unpack sol = integrator

  @inbounds if integrator.tdir * t < integrator.tdir * sol.t[1]
    if deriv == 0 && idxs === nothing
      return f.h(val, p, t)
    elseif idxs === nothing
      return f.h(val, p, t, Val{deriv})
    elseif deriv == 0
      return f.h(val, p, t; idxs = idxs)
    else
      return f.h(val, p, t, Val{deriv}; idxs = idxs)
    end
  elseif integrator.tdir * t <= integrator.tdir * sol.t[end] # Put equals back
    return sol.interp(val, t, idxs, Val{deriv}, f.integrator.p)
  end

  # history function is evaluated at a time point past the final time point of
  # the current solution
  f.isout = true

  # handle extrapolations at initial time point
  if integrator.t == sol.prob.tspan[1]
    return constant_extrapolant!(val, t, integrator, idxs, Val{deriv})
  else
    return OrdinaryDiffEq.current_interpolant!(val, t, f.integrator, idxs, Val{deriv})
  end
end
