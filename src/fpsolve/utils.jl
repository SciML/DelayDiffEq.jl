# construct solver for fixed-point iterations
function fpsolver(alg, u, uEltypeNoUnits, uBottomEltypeNoUnits, ::Val{iip}) where iip
  # no fixed-point iterations if the algorithm is constrained
  isconstrained(alg) && return

  uTolType = real(uBottomEltypeNoUnits)

  if iip
    du = similar(u)

    # could use integrator.cache.atmp if guaranteed to exist
    atmp = similar(u, uEltypeNoUnits)
  end

  # create cache
  if alg.fpsolve isa NLFunctional
    fpcache = iip ? FPFunctionalCache(du, atmp) : FPFunctionalConstantCache()
  elseif alg.fpsolve isa NLAnderson
    max_history = min(alg.fpsolve.max_history, alg.fpsolve.max_iter, length(u))

    Q = Matrix{uEltypeNoUnits}(undef, length(u), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    if iip
      Δus = [zero(u) for i in 1:max_history]
      duold = zero(u)
      uold = zero(u)
      fpcache = FPAndersonCache(du, duold, uold, atmp, Δus, Q, R, γs)
    else
      Δus = Vector{typeof(u)}(undef, max_history)
      fpcache = FPAndersonConstantCache(Δus, Q, R, γs)
    end
  end

  # create solver
  FPSolver{typeof(alg.fpsolve),iip,uTolType,typeof(fpcache)}(alg.fpsolve, one(uTolType), 10000, Convergence, fpcache)
end

# resize solver for fixed-point iterations
function fpsolve_resize!(integrator::DDEIntegrator, i)
  fpsolver = integrator.fpsolver

  if fpsolver !== nothing
    resize!(fpsolver, i)
  end

  nothing
end

function Base.resize!(fpsolver::FPSolver, i)
  if fpsolver.alg isa NLAnderson
    resize!(fpsolver.cache, fpsolver.alg, i)
  else
    resize!(fpsolver.cache, i)
  end

  nothing
end

function Base.resize!(fpcache::FPFunctionalCache, i::Int)
  resize!(fpcache.du, i)
  resize!(fpcache.atmp, i)

  nothing
end

function Base.resize!(fpcache::FPAndersonCache, fpalg::NLAnderson, i::Int)
  @unpack du, duold, uold, atmp, γs, Δus = fpcache

  resize!(du, i)
  resize!(duold, i)
  resize!(uold, i)
  resize!(atmp, i)

  # determine new maximum history
  max_history_old = length(Δus)
  max_history = min(fpalg.max_history, fpalg.max_iter, i)

  resize!(γs, max_history)
  resize!(Δus, max_history)
  if max_history > max_history_old
    for i in (max_history_old + 1):max_history
      Δus[i] = zero(uold)
    end
  end

  if max_history != max_history_old
    fpcache.Q = typeof(fpcache.Q)(undef, i, max_history)
    fpcache.R = typeof(fpcache.R)(undef, max_history, max_history)
  end

  nothing
end

function Base.resize!(fpcache::FPAndersonConstantCache, fpalg::NLAnderson, i::Int)
  @unpack γs, Δus = fpcache

  # determine new maximum history
  max_history_old = length(Δus)
  max_history = min(fpalg.max_history, fpalg.max_iter, i)

  resize!(γs, max_history)
  resize!(Δus, max_history)

  if max_history != max_history_old
    fpcache.Q = typeof(fpcache.Q)(undef, i, max_history)
    fpcache.R = typeof(fpcache.R)(undef, max_history, max_history)
  end

  nothing
end