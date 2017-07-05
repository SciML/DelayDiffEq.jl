# nlsolve optimizes subtypes of AbstractDifferentiableMultivariateFunction
mutable struct IterationFunction{I<:AbstractDDEIntegrator} <: AbstractDifferentiableMultivariateFunction
  integrator::I
  f!::IterationFunction{I} # Anderson acceleration needs field f!

  function IterationFunction{I}(integrator::I) where I <:AbstractDDEIntegrator
    f = new(integrator)
    f.f! = f
  end
end

function (f::IterationFunction)(x,fvec)
  integrator = f.integrator  

  # update u
  if typeof(integrator.u) <: AbstractArray
    copy!(integrator.u,x)
  else
    integrator.u = x[1] # x is always a vector
  end

  perform_step!(integrator,integrator.cache)

  # update output vector of residuals
  # in contrast to picard iteration we can not use integrator.resid
  # as output vector for different fixed-point iterations
  @. fvec = (integrator.u - x)/
    @muladd(integrator.fixedpoint_abstol + max(abs(integrator.u),abs(x))*integrator.fixedpoint_reltol)
 
  # special updates of ODE integrator after the first iteration necessary
  if integrator.first_iteration
    integrator.integrator.tprev = integrator.t
    integrator.integrator.t = integrator.t+integrator.dt
    integrator.integrator.uprev = integrator.u
    integrator.first_iteration = false
  end

  # update ODE integrator with values of DDE integrator
  integrator.integrator.u = integrator.u
  integrator.integrator.k = integrator.k
end
