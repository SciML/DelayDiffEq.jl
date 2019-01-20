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

function (f::HistoryFunction)(p, t, ::Type{Val{deriv}}=Val{0}; idxs=nothing) where deriv
    integrator = f.integrator

    @inbounds if integrator.tdir * t < integrator.tdir * f.sol.t[1]
        if deriv == 0 && idxs === nothing
            return f.h(p, t)
        elseif idxs === nothing
            return f.h(p, t, Val{deriv})
        elseif deriv == 0
            return f.h(p, t; idxs = idxs)
        else
            return f.h(p, t, Val{deriv}; idxs = idxs)
        end
    elseif integrator.tdir * t <= integrator.tdir * f.sol.t[end] # Put equals back
        return f.sol.interp(t, idxs, Val{deriv}, f.integrator.p)
    end

    # set boolean that indicates that history function was evaluated at time point past the
    # final time point of the current solution
    # NOTE: does not interfere with usual use of isout since this integrator is only used
    # for inter- and extrapolation of future values and saving of the solution but does not
    # affect whether time steps are accepted
    integrator.isout = true

    # handle extrapolations at initial time point
    if integrator.t == integrator.sol.prob.tspan[1]
        return constant_extrapolant(t, integrator, idxs, Val{deriv})
    else
        return OrdinaryDiffEq.current_interpolant(t, integrator, idxs, Val{deriv})
    end
end

function (f::HistoryFunction)(val, p, t, ::Type{Val{deriv}}=Val{0}; idxs=nothing) where deriv
    integrator = f.integrator

    @inbounds if integrator.tdir * t < integrator.tdir * f.sol.t[1]
        if deriv == 0 && idxs === nothing
            return f.h(val, p, t)
        elseif idxs === nothing
            return f.h(val, p, t, Val{deriv})
        elseif deriv == 0
            return f.h(val, p, t; idxs = idxs)
        else
            return f.h(val, p, t, Val{deriv}; idxs = idxs)
        end
    elseif integrator.tdir * t <= integrator.tdir * f.sol.t[end] # Put equals back
        return f.sol.interp(val, t, idxs, Val{deriv}, f.integrator.p)
    end

    # set boolean that indicates that history function was evaluated at time point past the
    # final time point of the current solution
    # NOTE: does not interfere with usual use of isout since this integrator is only used
    # for inter- and extrapolation of future values and saving of the solution but does not
    # affect whether time steps are accepted
    integrator.isout = true

    # handle extrapolations at initial time point
    if integrator.t == integrator.sol.prob.tspan[1]
        return constant_extrapolant!(val, t, integrator, idxs, Val{deriv})
    else
        return OrdinaryDiffEq.current_interpolant!(val, t, f.integrator, idxs, Val{deriv})
    end
end
