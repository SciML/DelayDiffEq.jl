function apply_callback!(integrator::DDEIntegrator,callback,cb_time=0,prev_sign=1)
  if cb_time != 0
    OrdinaryDiffEq.change_t_via_interpolation!(integrator,integrator.tprev+cb_time)
  end

  if callback.save_positions[1]
    savevalues!(integrator)
  end

  integrator.u_modified = true
  if prev_sign < 0 && !(typeof(callback.affect!) <: Void)
    callback.affect!(integrator)
  elseif !(typeof(callback.affect_neg!) <: Void)
    callback.affect_neg!(integrator)
  end
  if integrator.u_modified
    OrdinaryDiffEq.reeval_internals_due_to_modification!(integrator)
    push!(integrator.opts.d_discontinuities,compute_discontinuity_tree(integrator.sol.prob.lags,integrator.alg,integrator.t)...)
  end
  if callback.save_positions[2]
    savevalues!(integrator)
  end

end
