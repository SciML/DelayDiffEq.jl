"""
    HistoryFunction(h, integrator)

Wrap history function `h` and integrator `integrator` to create a common
interface for retrieving values at any time point with varying accuracy.

Before the initial time point of the dense solution values are calculated by history function
`h`, for time points in the time span of the solution the solution is interpolated, and
after the final time point of the solution an inter- or extrapolation of the current state
of integrator `integrator` is retrieved.
"""
mutable struct HistoryFunction{H,I<:AbstractODEIntegrator} <: Function
  h::H
  integrator::I
  isout::Bool
end

HistoryFunction(h, integrator::AbstractODEIntegrator) =
  HistoryFunction{typeof(h),typeof(integrator)}(h, integrator, false)

function (f::HistoryFunction)(p, t, ::Type{Val{deriv}}=Val{0}; idxs=nothing) where deriv
  @unpack h, integrator = f
  @unpack sol, tdir = integrator

    @inbounds if tdir * t < tdir * sol.t[1]
        if deriv == 0 && idxs === nothing
            return h(p, t)
        elseif idxs === nothing
            return h(p, t, Val{deriv})
        elseif deriv == 0
            return h(p, t; idxs = idxs)
        else
            return h(p, t, Val{deriv}; idxs = idxs)
        end
    elseif tdir * t <= tdir * sol.t[end] # Put equals back
        return sol.interp(t, idxs, Val{deriv}, p)
    end

    # set boolean that indicates that history function was evaluated at time point past the
    # final time point of the current solution
    f.isout = true

    # handle extrapolations at initial time point
    if integrator.t == sol.prob.tspan[1]
        return constant_extrapolant(t, integrator, idxs, Val{deriv})
    else
        return OrdinaryDiffEq.current_interpolant(t, integrator, idxs, Val{deriv})
    end
end

function (f::HistoryFunction)(val, p, t, ::Type{Val{deriv}}=Val{0}; idxs=nothing) where deriv
  @unpack h, integrator = f
  @unpack sol, tdir = integrator

  @inbounds if tdir * t < tdir * sol.t[1]
    if deriv == 0 && idxs === nothing
      return h(val, p, t)
    elseif idxs === nothing
      return h(val, p, t, Val{deriv})
    elseif deriv == 0
      return h(val, p, t; idxs = idxs)
    else
      return h(val, p, t, Val{deriv}; idxs = idxs)
    end
  elseif tdir * t <= tdir * sol.t[end] # Put equals back
    return sol.interp(val, t, idxs, Val{deriv}, p)
  end

  # set boolean that indicates that history function was evaluated at time point past the
  # final time point of the current solution
  f.isout = true

  # handle extrapolations at initial time point
  if integrator.t == sol.prob.tspan[1]
    return constant_extrapolant!(val, t, integrator, idxs, Val{deriv})
  else
    return OrdinaryDiffEq.current_interpolant!(val, t, integrator, idxs, Val{deriv})
  end
end
