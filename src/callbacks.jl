@inline function handle_callback_modifiers!(integrator::DDEIntegrator)
  integrator.reeval_fsal = true
  push!(integrator.opts.d_discontinuities,compute_discontinuity_tree(integrator.sol.prob.lags,integrator.alg,integrator.t)...)
end
