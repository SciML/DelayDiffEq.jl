# update integrator when u is modified by callbacks
@inline function handle_callback_modifiers!(integrator::DDEIntegrator)
    integrator.reeval_fsal = true # recalculate fsalfirst after applying step

    # update heap of discontinuities

    if typeof(integrator.sol.prob) <: ConstantLagDDEProblem
        #warn("ConstantLagDDEProblem is deprecated. Use DDEProblem instead.")
        neutral = false
        constant_lags = integrator.sol.prob.lags
    else
        neutral = integrator.sol.prob.neutral
        constant_lags = integrator.sol.prob.constant_lags
    end

    push!(integrator.opts.d_discontinuities,
          compute_discontinuity_tree(constant_lags, integrator.alg,
                                     integrator.t,integrator.sol.prob.tspan[2],neutral)...)
end
