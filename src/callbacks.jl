# update integrator when u is modified by callbacks
function handle_callback_modifiers!(integrator::DDEIntegrator)
    integrator.reeval_fsal = true # recalculate fsalfirst after applying step

    # update heap of discontinuities and time stops

    if typeof(integrator.sol.prob) <: ConstantLagDDEProblem
        #warn("ConstantLagDDEProblem is deprecated. Use DDEProblem instead.")
        neutral = false
        constant_lags = integrator.sol.prob.lags
    else
        neutral = integrator.sol.prob.neutral
        constant_lags = integrator.sol.prob.constant_lags
    end

    discontinuity_tree = compute_discontinuity_tree(constant_lags, integrator.alg,
                                                    integrator.t,
                                                    integrator.sol.prob.tspan[2], neutral)
    push!(integrator.opts.d_discontinuities, discontinuity_tree...)
    push!(integrator.opts.tstops, discontinuity_tree...)
end

"""
    reeval_internals_due_to_modification!(integrator::DDEIntegrator)

Recalculate interpolation data and update ODE integrator after changes by callbacks.
"""
function reeval_internals_due_to_modification!(integrator::DDEIntegrator)
    # update interpolation data of DDE integrator using old interpolation data
    # of ODE integrator in evaluation of history function that was calculated in
    # `perform_step!`
    OrdinaryDiffEq.ode_addsteps!(integrator, integrator.f, Val{true}, Val{true}, Val{true})

    # copy interpolation data to ODE integrator
    recursivecopy!(integrator.integrator.k, integrator.k)

    # move ODE integrator to new time interval of DDE integrator
    integrator.integrator.t = integrator.t
    integrator.integrator.dt = integrator.dt
    if typeof(integrator.integrator.u) <: AbstractArray
        recursivecopy!(integrator.integrator.u, integrator.u)
    else
        integrator.integrator.u = integrator.u
    end

    integrator.u_modified = false
end
