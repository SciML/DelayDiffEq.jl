function fpsolve!(integrator::DDEIntegrator, fpsolver::FPSolver{<:FPFunctionalConstantCache})
  @unpack f, t, p, k, uprev, dt, alg, cache = integrator
  @unpack κ, max_iter = fpsolver
  ode_integrator = integrator.integrator

  # ODE integrator caches state u at next time point
  # DDE integrator contains updated state u₊ at next time point

  # precalculations
  η = fpsolver.ηold

  # update ODE integrator to next time interval together with correct interpolation
  advance_ode_integrator!(integrator)

  # fixed point iteration
  local ndu
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1

    # calculate next step
    OrdinaryDiffEq.perform_step!(integrator, integrator.cache, true)

    # compute norm of residuals
    iter > 1 && (nduprev = ndu)
    atmp = OrdinaryDiffEq.calculate_residuals(ode_integrator.u, integrator.u,
                                              integrator.opts.abstol,
                                              integrator.opts.reltol,
                                              integrator.opts.internalnorm, t)
    ndu = integrator.opts.internalnorm(atmp, t)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndu / nduprev
      ( diverge = θ > 1 ) && ( fpsolver.status = Divergence )
      ( veryslowconvergence = ndu * θ^(max_iter - iter) > κ * (1 - θ) ) && ( fpsolver.status = VerySlowConvergence )
      if diverge || veryslowconvergence
        break
      end
    end

    # complete interpolation data of DDE integrator for time interval [t, t+dt]
    # and copy it to ODE integrator
    # has to be done before updates to ODE integrator, otherwise history function
    # is incorrect
    if iscomposite(alg)
      addsteps!(k, t, uprev, integrator.u, dt, f, p, cache.caches[integrator.cache.current],
                false, true, true)
    else
      addsteps!(k, t, uprev, integrator.u, dt, f, p, cache, false, true, true)
    end
    @inbounds for i in 1:length(k)
      copyat_or_push!(ode_integrator.k, i, k[i])
    end

    # update state of the dummy ODE solver
    ode_integrator.u = integrator.u

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndu < κ && (iter > 1 || iszero(ndu))
      # fixed-point iteration converges
      fpsolver.status = η < fpsolver.fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end
  end

  integrator.force_stepfail = fail_convergence || integrator.force_stepfail
  fpsolver.ηold = η
  fpsolver.fp_iters = iter

  nothing
end

function fpsolve!(integrator::DDEIntegrator, fpsolver::FPSolver{<:FPFunctionalCache})
  @unpack f, t, p, k, uprev, dt, alg = integrator
  @unpack κ, max_iter, cache = fpsolver
  @unpack resid = cache
  ode_integrator = integrator.integrator

  # ODE integrator caches state u at next time point
  # DDE integrator contains updated state u₊ at next time point

  # precalculations
  η = fpsolver.ηold

  # update ODE integrator to next time interval together with correct interpolation
  advance_ode_integrator!(integrator)

  # fixed-point iteration without Newton
  local ndu
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1

    # calculate next step
    OrdinaryDiffEq.perform_step!(integrator, integrator.cache, true)

    # compute norm of residuals
    iter > 1 && (nduprev = ndu)
    OrdinaryDiffEq.calculate_residuals!(resid, ode_integrator.u, integrator.u,
                                        integrator.opts.abstol, integrator.opts.reltol,
                                        integrator.opts.internalnorm, t)
    ndu = integrator.opts.internalnorm(resid, t)

    # check divergence (not in initial step)
    if iter > 1
      θ = ndu / nduprev
      ( diverge = θ > 1 ) && ( fpsolver.status = Divergence )
      ( veryslowconvergence = ndu * θ^(max_iter - iter) > κ * (1 - θ) ) && ( fpsolver.status = VerySlowConvergence )
      if diverge || veryslowconvergence
        break
      end
    end

    # complete interpolation data of DDE integrator for time interval [t, t+dt]
    # and copy it to ODE integrator
    # has to be done before updates to ODE integrator, otherwise history function
    # is incorrect
    if iscomposite(alg)
      addsteps!(k, t, uprev, integrator.u, dt, f, p, cache.caches[integrator.cache.current],
                false, true, true)
    else
      addsteps!(k, t, uprev, integrator.u, dt, f, p, cache, false, true, true)
    end
    @inbounds for i in 1:length(k)
      copyat_or_push!(ode_integrator.k, i, k[i])
    end

    # update state of the dummy ODE solver
    recursivecopy!(ode_integrator.u, integrator.u)

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndu < κ && (iter > 1 || iszero(ndu))
      # fixed-point iteration converges
      fpsolver.status = η < fpsolver.fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end
  end

  integrator.force_stepfail = fail_convergence || integrator.force_stepfail
  fpsolver.ηold = η
  fpsolver.fp_iters = iter

  nothing
end
