## preamble!

function DiffEqBase.initialize_cache!(fpcache::Union{FPAndersonConstantCache,FPAndersonCache},
                                      fpsolver::FPSolver{<:NLAnderson}, integrator)
  fpcache.history = 0
  
  nothing
end

## apply_step!

function DiffEqBase.apply_step!(fpsolver::FPSolver{<:NLAnderson,false}, integrator)
  @unpack iter,alg,cache = fpsolver
  @unpack aa_start = alg

  # perform Anderson acceleration
  if iter == aa_start
    # update cached values for next step of Anderson acceleration
    cache.gzprev = integrator.u
    cache.dzprev = @. integrator.u - integrator.integrator.u
  elseif iter > aa_start
    # actually update the next iterate
    integrator.u = anderson(integrator.u, cache, integrator)

    # update the initial stages
    addsteps!(integrator.k, integrator.t, integrator.uprev, integrator.u, integrator.dt,
              integrator.f, integrator.p, integrator.cache, true)
  end

  # apply the step, i.e., update the interpolation
  _apply_step!(fpsolver, integrator)

  nothing
end

function DiffEqBase.apply_step!(fpsolver::FPSolver{<:NLAnderson,true}, integrator)
  @unpack iter,alg,cache = fpsolver
  @unpack aa_start = alg

  # perform Anderson acceleration
  if iter == aa_start
    # update cached values for next step of Anderson acceleration
    @.. cache.gzprev = integrator.u
    @.. cache.dzprev = integrator.u - integrator.integrator.u
  elseif iter > aa_start
    # actually update the next stage
    anderson!(integrator.u, fpsolver.cache, integrator)

    # update the initial stages
    addsteps!(integrator.k, integrator.t, integrator.uprev, integrator.u, integrator.dt,
              integrator.f, integrator.p, integrator.cache, true)
  end

  # apply the step, i.e., update the interpolation
  _apply_step!(fpsolver, integrator)

  nothing
end

## perform_step!

DiffEqBase.perform_step!(fpsolver::FPSolver{<:NLAnderson}, integrator) =
  OrdinaryDiffEq.perform_step!(integrator, integrator.cache, true)

## norm_of_residuals

function DiffEqBase.norm_of_residuals(fpsolver::FPSolver{<:NLAnderson,false}, integrator)
  @unpack t,opts = integrator
  @unpack cache = fpsolver

  zprev = integrator.integrator.u
  z = integrator.u
  cache.dz = @. z - zprev
  atmp = calculate_residuals(cache.dz, zprev, z, opts.abstol, opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end

function DiffEqBase.norm_of_residuals(fpsolver::FPSolver{<:NLAnderson,true}, integrator)
  @unpack t,opts = integrator
  @unpack dz,atmp = fpsolver.cache

  zprev = integrator.integrator.u
  z = integrator.u
  @.. dz = z - zprev
  calculate_residuals!(atmp, dz, zprev, z, opts.abstol, opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end

## resize!

function Base.resize!(fpcache::Union{FPAndersonCache,FPAndersonConstantCache},
                      fpsolver::FPSolver{<:NLAnderson}, integrator, i::Int)
  resize!(fpcache, fpsolver.alg, i)
end

function Base.resize!(fpcache::FPAndersonCache, fpalg::NLAnderson, i::Int)
  @unpack gzprev, Δgzs = fpcache
  
  resize!(fpcache.dz, i)
  resize!(fpcache.atmp, i)
  resize!(gzprev, i)
  resize!(fpcache.dzprev, i)

  # update history of Anderson cache
  max_history_old = length(Δgzs)
  DiffEqBase.resize_anderson_history!(fpcache, fpalg, i)

  max_history = length(Δgzs)
  if max_history > max_history_old
    for i in (max_history_old + 1):max_history
      Δgzs[i] = zero(gzprev)
    end
  end

  if max_history != max_history_old
    nlcache.Q = typeof(nlcache.Q)(undef, i, max_history)
    nlcache.R = typeof(nlcache.R)(undef, max_history, max_history)
  end

  nothing
end

Base.resize!(fpcache::FPAndersonConstantCache, fpalg::NLAnderson, i::Int) =
  DiffEqBase.resize_anderson_history!(fpcache, fpalg, i)