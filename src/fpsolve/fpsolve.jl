## default implementations

function DiffEqBase.norm_of_residuals(fpsolver::AbstractFPSolver, integrator)
  @unpack t,opts = integrator

  atmp = calculate_residuals(integrator.integrator.u, integrator.u, opts.abstol,
                             opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end

function DiffEqBase.postamble!(fpsolver::AbstractFPSolver, integrator)
  fail_convergence = nlsolvefail(fpsolver)
  # TODO: update statistics
  integrator.force_stepfail = fail_convergence || integrator.force_stepfail

  # we do not have any convergence guarantees on integrator.u, hence
  # we reset it to integrator.integrator.u and the corresponding interpolation
  reset_integrator!(integrator)
end
