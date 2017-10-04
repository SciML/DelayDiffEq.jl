"""
    HistoryFunction(h, sol, integrator)

Wrap history function `h`, solution `sol`, and integrator `integrator` to create a common
interface for retrieving values at any time point with varying accuracy.

Before the initial time point of solution `sol` values are calculated by history function
`h`, for time points in the time span of `sol` interpolated values of `sol` are returned,
and after the final time point of `sol` an inter- or extrapolation of the current state
of integrator `integrator` is retrieved.
"""
struct HistoryFunction{F1,F2,F3<:ODEIntegrator} <: Function
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
    # NOTE: does not interfere with usual use of isout since this integrator is only used
    # for inter- and extrapolation of future values and saving of the solution but does not
    # affect whether time steps are accepted
    integrator.isout = true

    # handle extrapolations at initial time point
    f.sol.t[end] == integrator.sol.prob.tspan[1] && return constant_extrapolant(t, integrator, idxs, deriv)

    return OrdinaryDiffEq.current_interpolant(t, integrator, idxs, deriv)
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
    # NOTE: does not interfere with usual use of isout since this integrator is only used
    # for inter- and extrapolation of future values and saving of the solution but does not
    # affect whether time steps are accepted
    integrator.isout = true

    # handle extrapolations at initial time point
    f.sol.t[end] == integrator.sol.prob.tspan[1] && return constant_extrapolant!(val, t, integrator, idxs, deriv)

    return OrdinaryDiffEq.current_interpolant!(val, t, f.integrator, idxs, deriv)
end
