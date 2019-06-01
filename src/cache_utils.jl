function DiffEqBase.unwrap_cache(integrator::DDEIntegrator, is_stiff)
  alg   = integrator.alg
  cache = integrator.cache
  iscomp = alg isa CompositeAlgorithm
  if !iscomp
    return cache
  elseif alg.choice_function isa AutoSwitch
    num = is_stiff ? 2 : 1
    return cache.caches[num]
  else
    return cache.caches[integrator.cache.current]
  end
end
