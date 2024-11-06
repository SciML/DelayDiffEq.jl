OrdinaryDiffEqNonlinearSolve.initial_η(fpsolver::FPSolver, integrator::DDEIntegrator) = fpsolver.ηold

OrdinaryDiffEqNonlinearSolve.apply_step!(fpsolver::FPSolver, integrator::DDEIntegrator) = nothing

function DiffEqBase.postamble!(fpsolver::FPSolver, integrator::DDEIntegrator)
    integrator.stats.nfpiter += fpsolver.iter

    if OrdinaryDiffEqNonlinearSolve.nlsolvefail(fpsolver)
        integrator.stats.nfpconvfail += 1
    end
    integrator.force_stepfail = OrdinaryDiffEqNonlinearSolve.nlsolvefail(fpsolver) ||
                                integrator.force_stepfail

    nothing
end
