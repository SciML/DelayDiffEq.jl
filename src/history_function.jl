"""
    HistoryFunction(h, sol, integrator)

Wrap history function `h`, solution `sol`, and integrator `integrator` to create a common
interface for retrieving values at any time point with varying accuracy.

Before the initial time point of solution `sol` values are calculated by history function
`h`, for time points in the time span of `sol` interpolated values of `sol` are returned,
and after the final time point of `sol` an inter- or extrapolation of the current state
of integrator `integrator` is retrieved.
"""
struct HistoryFunction{F1,F2,F3} <: Function
    h::F1
    sol::F2
    integrator::F3
end

function (f::HistoryFunction)(t, deriv::Type=Val{0}, idxs=nothing)
    @inbounds if t < f.sol.t[1]
        if typeof(idxs) <: Void
            return f.h(t)
        else
            return f.h(t, idxs)
        end
    elseif t <= f.sol.t[end] # Put equals back
        return f.sol.interp(t, idxs, deriv)
    end

    integrator = f.integrator

    # set boolean that indicates that history function was evaluated at time point past the
    # final time point of the current solution
    integrator.isout = true

    # handle extrapolations in initial time step (otherwise dt = 0 should not occur)
    integrator.dt == 0 && return constant_extrapolation(t, integrator, idxs, deriv)

    if typeof(integrator.alg) <: Union{Vern6,Vern7,Vern8,Vern9} &&
        length(integrator.k) == integrator.kshortsize
        # use Hermite interpolant as fallback to avoid implicit calculation
        return hermite_interpolant(t, integrator, idxs, deriv)
    else
        return OrdinaryDiffEq.current_interpolant(t, integrator, idxs, deriv)
    end
end

function (f::HistoryFunction)(val, t, deriv::Type=Val{0}, idxs=nothing)
    @inbounds if t < f.sol.t[1]
        if typeof(idxs) <: Void
            return f.h(val, t)
        else
            return f.h(val, t, idxs)
        end
    elseif t <= f.sol.t[end] # Put equals back
        return f.sol.interp(val, t, idxs, deriv)
    end

    integrator = f.integrator

    # set boolean that indicates that history function was evaluated at time point past the
    # final time point of the current solution
    integrator.isout = true

    # handle extrapolations in initial time step (otherwise dt = 0 should not occur)
    integrator.dt == 0 && return constant_extrapolation!(val, t, integrator, idxs, deriv)

    if typeof(integrator.alg) <: Union{Vern6,Vern7,Vern8,Vern9} &&
        length(integrator.k) == integrator.kshortsize
        # use Hermite interpolant as fallback to avoid implicit calculation
        return hermite_interpolant!(val, t, integrator, idxs, deriv)
    else
        return OrdinaryDiffEq.current_interpolant!(val, t, f.integrator, idxs, deriv)
    end
end
