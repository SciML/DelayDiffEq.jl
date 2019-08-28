function get_fpsolver(integrator::DEIntegrator)
  isdefined(integrator.cache, :nlsolver) || return
    
  integrator.cache.nlsolver
end

## resize

function resize_fpsolver!(integrator::DEIntegrator, i)
  fpsolver = get_fpsolver(integrator)
  
  fpsolver === nothing && return

  resize!(fpsolver, integrator, i)
  
  nothing
end
  
Base.resize!(fpsolver::AbstractFPSolver, integrator, i::Int) =
  resize!(DiffEqBase.get_cache(fpsolver), fpsolver, integrator, i)