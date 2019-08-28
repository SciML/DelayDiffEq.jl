## build_nlsolver

function build_fpsolver(alg, fpalg::Union{NLFunctional,NLAnderson}, u, uEltypeNoUnits,
                        uBottomEltypeNoUnits, ::Val{true})
  # no fixed-point iterations if the algorithm is constrained
  isconstrained(alg) && return

  # define unitless types
  uTolType = real(uBottomEltypeNoUnits)

  # build cache for fixed-point solver
  atmp = similar(u, uEltypeNoUnits)

  if fpalg isa NLFunctional
    cache = FPFunctionalCache(atmp)
  elseif fpalg isa NLAnderson
    dz = similar(u)
    dzprev = similar(u)
    gzprev = similar(u)

    max_history = min(fpalg.max_history, fpalg.maxiters, length(u))
    Δgzs = [zero(u) for i in 1:max_history]
    Q = Matrix{uEltypeNoUnits}(undef, length(u), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    droptol = fpalg.droptol === nothing ? nothing : uEltypeNoUnits(fpalg.droptol)

    cache = FPAndersonCache(dz, atmp, gzprev, dzprev, Δgzs, Q, R, γs, 0, droptol)
  end

  # build non-linear solver
  η = one(uTolType)
  ndu = one(uTolType)

  FPSolver{typeof(fpalg),true,uTolType,typeof(cache)}(
    fpalg, uTolType(fpalg.κ), η, ndu, uTolType(fpalg.fast_convergence_cutoff),
    fpalg.maxiters, 10_000, Convergence, cache)  
end

function build_fpsolver(alg, fpalg::Union{NLFunctional,NLAnderson}, u, uEltypeNoUnits,
                        uBottomEltypeNoUnits, ::Val{false})
  # no fixed-point iterations if the algorithm is constrained
  isconstrained(alg) && return

  # define unitless types
  uTolType = real(uBottomEltypeNoUnits)

  # create cache of fixed-point solver
  if fpalg isa NLFunctional
    cache = FPFunctionalConstantCache()
  elseif fpalg isa NLAnderson
    dz = u
    dzprev = u
    gzprev = u
    
    max_history = min(fpalg.max_history, fpalg.maxiters, length(z))
    Δgzs = Vector{typeof(u)}(undef, max_history)
    Q = Matrix{uEltypeNoUnits}(undef, length(u), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    droptol = fpalg.droptol === nothing ? nothing : uEltypeNoUnits(fpalg.droptol)

    cache = FPAndersonConstantCache(dz, gzprev, dzprev, Δgzs, Q, R, γs, 0, droptol)
  end

  # build non-linear solver
  η = one(uTolType)
  ndu = one(uTolType)
  
  FPSolver{typeof(fpalg),false,uTolType,typeof(cache)}(
    fpalg, uTolType(fpalg.κ), η, ndu, uTolType(fpalg.fast_convergence_cutoff),
    fpalg.maxiters, 10_000, Convergence, cache)
end

## _apply_step!

function DiffEqBase._apply_step!(fpsolver::FPSolver, integrator)
  # update ODE integrator to next time interval together with correct interpolation
  if fpsolver.iter > 0
    update_ode_integrator!(integrator)
  end

  # reset indicator for failed non-linear solvers
  integrator.force_stepfail = false

  # update statistics
  fpsolver.iter += 1
  # TODO: update statistics

  nothing
end