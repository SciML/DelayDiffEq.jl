## perform_step!

DiffEqBase.perform_step!(::FPSolver{<:NLFunctional}, integrator) =
  OrdinaryDiffEq.perform_step!(integrator, integrator.cache, true)

## norm_of_residuals

function DiffEqBase.norm_of_residuals(fpsolver::FPSolver{<:NLFunctional,true},
                                      integrator)
  @unpack t,opts = integrator
  @unpack atmp = fpsolver.cache

  calculate_residuals!(atmp, integrator.integrator.u, integrator.u,
                       opts.abstol, opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end

## resize!

function Base.resize!(fpcache::FPFunctionalCache, i::Int)
  resize!(fpcache.atmp, i)

  nothing
end