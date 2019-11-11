OrdinaryDiffEq.initial_η(fpsolver::FPSolver, integrator::DDEIntegrator) = fpsolver.ηold

OrdinaryDiffEq.apply_step!(fpsolver::FPSolver, integrator::DDEIntegrator) = nothing

function DiffEqBase.postamble!(fpsolver::FPSolver, integrator::DDEIntegrator)
  integrator.destats.nfpiter += fpsolver.iter

  if OrdinaryDiffEq.nlsolvefail(fpsolver)
    integrator.destats.nfpconvfail += 1
  end
  integrator.force_stepfail = OrdinaryDiffEq.nlsolvefail(fpsolver) || integrator.force_stepfail

  nothing
end
