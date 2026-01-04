function OrdinaryDiffEqNonlinearSolve.initial_η(fpsolver::FPSolver, integrator::DDEIntegrator)
    return fpsolver.ηold
end

OrdinaryDiffEqCore.apply_step!(fpsolver::FPSolver, integrator::DDEIntegrator) = nothing

function SciMLBase.postamble!(fpsolver::FPSolver, integrator::DDEIntegrator)
    integrator.stats.nfpiter += fpsolver.iter

    if OrdinaryDiffEqNonlinearSolve.nlsolvefail(fpsolver)
        integrator.stats.nfpconvfail += 1
    end
    integrator.force_stepfail = OrdinaryDiffEqNonlinearSolve.nlsolvefail(fpsolver) ||
        integrator.force_stepfail

    return nothing
end
