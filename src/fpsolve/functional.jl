function fpsolve!(fpsolver::FPSolver{<:Union{NLFunctional,NLAnderson},false},
                  integrator::DDEIntegrator)
  @unpack f, t, p, k, uprev, dt, alg, cache = integrator
  ode_integrator = integrator.integrator

  @unpack κ, max_iter, fast_convergence_cutoff = fpsolver.alg
  fpcache = fpsolver.cache

  if fpcache isa FPAndersonConstantCache
    @unpack aa_start,droptol = fpsolver.alg
    @unpack Δus,Q,R,γs = fpcache
    local duold, uold
  end

  # ODE integrator caches state u at next time point
  # DDE integrator contains updated state u₊ at next time point

  # precalculations
  η = fpsolver.ηold
  if fpcache isa FPAndersonConstantCache
    history = 0
    max_history = length(Δus)
  end

  # fixed point iteration
  local ndu
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1

    # update ODE integrator to next time interval together with correct interpolation
    if iter == 1
      advance_ode_integrator!(integrator)
    else
      # for Anderson acceleration force recomputation of all interpolation data
      update_ode_integrator!(integrator,
                             fpsolver.alg isa NLAnderson && iter > aa_start + 1)
    end

    # calculate next step
    OrdinaryDiffEq.perform_step!(integrator, cache, true)

    # compute norm of residuals
    iter > 1 && (nduprev = ndu)
    du = integrator.u .- ode_integrator.u
    atmp = OrdinaryDiffEq.calculate_residuals(du, ode_integrator.u, integrator.u,
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

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndu < κ && (iter > 1 || iszero(ndu))
      # fixed-point iteration converges
      fpsolver.status = η < fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end

    # perform Anderson acceleration
    if fpcache isa FPAndersonConstantCache && iter < max_iter
      if iter == aa_start
        # update cached values for next step of Anderson acceleration
        duold = du
        uold = integrator.u
      elseif iter > aa_start
        # increase size of history
        history += 1

        # remove oldest history if maximum size is exceeded
        if history > max_history
          # circularly shift differences of u
          for i in 1:(max_history-1)
            Δus[i] = Δus[i + 1]
          end

          # delete left-most column of QR decomposition
          DiffEqBase.qrdelete!(Q, R, max_history)

          # update size of history
          history = max_history
        end

        # update history of differences of u
        Δus[history] = @. integrator.u - uold

        # replace/add difference of residuals as right-most column to QR decomposition
        DiffEqBase.qradd!(Q, R, DiffEqBase._vec(du .- duold), history)

        # update cached values
        duold = du
        uold = integrator.u

        # define current Q and R matrices
        Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

        # check condition (TODO: incremental estimation)
        if droptol !== nothing
          while cond(R) > droptol && history > 1
            DiffEqBase.qrdelete!(Q, R, history)
            history -= 1
            Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
          end
        end

        # solve least squares problem
        γscur = view(γs, 1:history)
        ldiv!(Rcur, mul!(γscur, Qcur', DiffEqBase._vec(du)))

        # update next iterate
        for i in 1:history
          integrator.u = @. integrator.u - γs[i] * Δus[i]
        end

        # update norm of residuals
        du = integrator.u .- uold .+ du
        atmp = OrdinaryDiffEq.calculate_residuals(du, ode_integrator.u, integrator.u,
                                                  integrator.opts.abstol,
                                                  integrator.opts.reltol,
                                                  integrator.opts.internalnorm, t)
        ndu = integrator.opts.internalnorm(atmp, t)
      end
    end
  end

  integrator.force_stepfail = fail_convergence || integrator.force_stepfail
  fpsolver.ηold = η
  fpsolver.fp_iters = iter

  nothing
end

function fpsolve!(fpsolver::FPSolver{<:Union{NLFunctional,NLAnderson},true}, integrator::DDEIntegrator)
  @unpack f, t, p, k, uprev, dt, alg, cache = integrator
  ode_integrator = integrator.integrator

  @unpack κ, max_iter, fast_convergence_cutoff = fpsolver.alg
  fpcache = fpsolver.cache
  @unpack du,atmp = fpcache

  if fpcache isa FPAndersonCache
    @unpack aa_start,droptol = fpsolver.alg
    @unpack duold,uold,Δus,Q,R,γs = fpcache
  end

  # ODE integrator caches state u at next time point
  # DDE integrator contains updated state u₊ at next time point

  # precalculations
  η = fpsolver.ηold
  if fpcache isa FPAndersonCache
    history = 0
    max_history = length(Δus)
  end

  # fixed-point iteration without Newton
  local ndu
  fail_convergence = true
  iter = 0
  while iter < max_iter
    iter += 1

    # update ODE integrator to next time interval together with correct interpolation
    if iter == 1
      advance_ode_integrator!(integrator)
    else
      # for Anderson acceleration force recomputation of all interpolation data
      update_ode_integrator!(integrator,
                             fpsolver.alg isa NLAnderson && iter > aa_start + 1)
    end

    # calculate next step
    OrdinaryDiffEq.perform_step!(integrator, cache, true)

    # compute norm of residuals
    iter > 1 && (nduprev = ndu)
    @.. du = integrator.u - ode_integrator.u
    OrdinaryDiffEq.calculate_residuals!(atmp, du, ode_integrator.u, integrator.u,
                                        integrator.opts.abstol, integrator.opts.reltol,
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

    # check stopping criterion
    iter > 1 && (η = θ / (1 - θ))
    if η * ndu < κ && (iter > 1 || iszero(ndu))
      # fixed-point iteration converges
      fpsolver.status = η < fast_convergence_cutoff ? FastConvergence : Convergence
      fail_convergence = false
      break
    end

    # perform Anderson acceleration
    if fpcache isa FPAndersonConstantCache && iter < max_iter
      if iter == aa_start
        # update cached values for next step of Anderson acceleration
        @.. duold = du
        @.. uold = integrator.u
      elseif iter > aa_start
        # increase size of history
        history += 1

        # remove oldest history if maximum size is exceeded
        if history > max_history
          # circularly shift differences of u
          ptr = Δus[1]
          for i in 1:(max_history-1)
            Δus[i] = Δus[i + 1]
          end
          Δus[max_history] = ptr

          # delete left-most column of QR decomposition
          DiffEqBase.qrdelete!(Q, R, max_history)

          # update size of history
          history = max_history
        end

        # update history of differences of u
        @.. Δus[history] = integrator.u - uold

        # replace/add difference of residuals as right-most column to QR decomposition
        @.. duold = du - duold
        DiffEqBase.qradd!(Q, R, vec(duold), history)

        # update cached values
        @.. duold = du
        @.. uold = integrator.u

        # define current Q and R matrices
        Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

        # check condition (TODO: incremental estimation)
        if droptol !== nothing
          while cond(R) > droptol && history > 1
            DiffEqBase.qrdelete!(Q, R, history)
            history -= 1
            Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
          end
        end

        # solve least squares problem
        γscur = view(γs, 1:history)
        ldiv!(Rcur, mul!(γscur, Qcur', vec(du)))

        # update next iterate
        for i in 1:history
          @.. integrator.u = integrator.u - γs[i] * Δus[i]
        end

        # update norm of residuals
        @.. du = integrator.u - uold + du
        OrdinaryDiffEq.calculate_residuals!(atmp, du, ode_integrator.u, integrator.u,
                                            integrator.opts.abstol,
                                            integrator.opts.reltol,
                                            integrator.opts.internalnorm, t)
        ndu = integrator.opts.internalnorm(atmp, t)
      end
    end
  end

  integrator.force_stepfail = fail_convergence || integrator.force_stepfail
  fpsolver.ηold = η
  fpsolver.fp_iters = iter

  nothing
end
