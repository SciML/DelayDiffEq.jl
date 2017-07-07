# update integrator when u is modified by callbacks
@inline function handle_callback_modifiers!(integrator::DDEIntegrator)
    integrator.reeval_fsal = true # recalculate fsalfirst after applying step

    # update heap of discontinuities
    push!(integrator.opts.d_discontinuities,
          compute_discontinuity_tree(integrator.sol.prob.lags, integrator.alg,
                                     integrator.t)...)
end
