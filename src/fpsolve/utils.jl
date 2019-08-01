# construct solver for fixed-point iterations
function fpsolver(alg, u, uEltypeNoUnits, uBottomEltypeNoUnits, ::Val{iip}) where iip
  # return dummy solver if the algorithm is constrained
  isconstrained(alg) && return FPNone()

  @unpack κ, fast_convergence_cutoff, max_iter = alg.fpsolve

  uTolType = real(uBottomEltypeNoUnits)

  # create cache
  if alg.fpsolve isa FPFunctional
    if iip
      z = similar(u, uEltypeNoUnits)
      fpcache = FPFunctionalCache(z)
    else
      fpcache = FPFunctionalConstantCache()
    end
  end

  # create solver
  FPSolver{typeof(fpcache),uTolType,typeof(κ),typeof(fast_convergence_cutoff)}(one(uTolType), κ, max_iter, 10000, Convergence, fast_convergence_cutoff, fpcache)
end

# perform fixed-point iteration
fpsolve!(integrator::DDEIntegrator) = fpsolve!(integrator, integrator.fpsolver)
fpsolve!(integrator, ::FPNone) = (@show integrator.t, integrator.tprev, integrator.dt, integrator.integrator.t, integrator.integrator.tprev, integrator.integrator.sol.t; error("Tried to perform fixed-point iteration although the algorithm is constrained. This is strictly forbidden."))

# resize solver for fixed-point iterations
fpsolve_resize!(integrator::DDEIntegrator, i) =
  fpsolve_resize!(integrator, integrator.fpsolver, i)

fpsolve_resize!(integrator, fpsolver, i) = nothing
function fpsolve_resize!(integrator, fpsolver::FPSolver{<:FPFunctionalCache}, i)
  resize!(fpsolver.cache.resid, i)
  nothing
end
