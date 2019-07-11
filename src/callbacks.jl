# update integrator when u is modified by callbacks
function handle_callback_modifiers!(integrator::DDEIntegrator)
    integrator.reeval_fsal = true # recalculate fsalfirst after applying step

    # update heap of discontinuities
    # discontinuity is assumed to be of order 0, i.e. solution x is discontinuous
    push!(integrator.opts.d_discontinuities, Discontinuity(integrator.t, 0))
end

"""
    reeval_internals_due_to_modification!(integrator::DDEIntegrator)

Recalculate interpolation data and update ODE integrator after changes by callbacks.
"""
function DiffEqBase.reeval_internals_due_to_modification!(integrator::DDEIntegrator,
            x::Type{Val{not_initialization}} = Val{true}) where not_initialization

    # update interpolation data of DDE integrator using old interpolation data
    # of ODE integrator in evaluation of history function that was calculated in
    # `perform_step!`
    if not_initialization
        addsteps!(integrator, integrator.f, true, true, true)
        # copy interpolation data to ODE integrator
        recursivecopy!(integrator.integrator.k, integrator.k)
    end

    # move ODE integrator to new time interval of DDE integrator
    integrator.integrator.t = integrator.t
    integrator.integrator.dt = integrator.dt
    if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.integrator.u, integrator.u)
    else
        integrator.integrator.u = integrator.u
    end

    integrator.u_modified = false
end
