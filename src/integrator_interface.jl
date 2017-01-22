function savevalues!(integrator::DDEIntegrator)
  integrator.integrator.u = integrator.u
  integrator.integrator.k = integrator.k
  integrator.integrator.t = integrator.t
  OrdinaryDiffEq.ode_addsteps!(integrator.integrator,integrator.f)
  savevalues!(integrator.integrator)
end

function postamble!(integrator::DDEIntegrator)
  integrator.integrator.u = integrator.u
  integrator.integrator.k = integrator.k
  integrator.integrator.t = integrator.t
  savevalues!(integrator.integrator)
end

function perform_step!(integrator::DDEIntegrator)
  integrator.integrator.uprev = integrator.uprev
  integrator.integrator.tprev = integrator.tprev
  integrator.integrator.fsalfirst = integrator.fsalfirst
  integrator.integrator.t = integrator.t
  integrator.integrator.dt = integrator.dt

  # if dt>max lag, then it's explicit so use Picard iteration
  if integrator.dt >minimum(integrator.prob.lags)

    # the is done to correct the extrapolation
    t_prev_cache = integrator.tprev
    t_cache = integrator.t
    uprev_cache = integrator.uprev

    numiters = 1
    while true
      if typeof(integrator.u) <: AbstractArray
        copy!(integrator.u_cache,integrator.u)
      else
        integrator.u_cache = integrator.u
      end
      perform_step!(integrator,integrator.cache)

      if typeof(integrator.resid) <: AbstractArray
        integrator.resid .= (integrator.u .- integrator.u_cache)./(integrator.picardabstol .+ max.(abs.(integrator.u),abs.(integrator.u_cache))*integrator.picardreltol)
      else
        integrator.resid = (integrator.u .- integrator.u_cache)./(integrator.picardabstol .+ max.(abs.(integrator.u),abs.(integrator.u_cache))*integrator.picardreltol)
      end

      picardEEst = integrator.picardnorm(integrator.resid)
      if picardEEst < 1 || numiters > integrator.max_picard_iters
        break
      end
      if numiters == 1
        integrator.integrator.tprev = integrator.t
        integrator.integrator.t = integrator.t+integrator.dt
        integrator.integrator.uprev = integrator.u
      end
      integrator.integrator.u = integrator.u
      integrator.integrator.k = integrator.k
      numiters += 1
    end
    integrator.t = t_cache
    integrator.integrator.t = t_cache
    integrator.tprev = t_prev_cache
    integrator.uprev = uprev_cache

  else # no iterations
    perform_step!(integrator,integrator.cache)
  end

  #integrator.u = integrator.integrator.u
  #integrator.fsallast = integrator.integrator.fsallast
  #if integrator.opts.adaptive
  #  integrator.EEst = integrator.integrator.EEst
  #end
end

function initialize!(dde_int::DDEIntegrator)
  initialize!(dde_int,dde_int.cache,dde_int.f)
  initialize!(dde_int.integrator,dde_int.cache,dde_int.f)
end

u_cache(integrator::DDEIntegrator) = u_cache(integrator.cache)
du_cache(integrator::DDEIntegrator)= du_cache(integrator.cache)
full_cache(integrator::DDEIntegrator) = chain(u_cache(integrator),du_cache(integrator.cache))

resize!(integrator::DDEIntegrator,i::Int) = resize!(integrator,integrator.cache,i)
function resize!(integrator::DDEIntegrator,cache,i)
  for c in full_cache(integrator)
    resize!(c,i)
  end
end

function resize!(integrator::DDEIntegrator,cache::Union{Rosenbrock23Cache,Rosenbrock32Cache},i)
  for c in full_cache(integrator)
    resize!(c,i)
  end
  for c in vecu_cache(integrator.cache)
    resize!(c,i)
  end
  Jvec = vec(cache.J)
  cache.J = reshape(resize!(Jvec,i*i),i,i)
  Wvec = vec(cache.W)
  cache.W = reshape(resize!(Wvec,i*i),i,i)
end

function resize!(integrator::DDEIntegrator,cache::Union{ImplicitEulerCache,TrapezoidCache},i)
  for c in full_cache(integrator)
    resize!(c,i)
  end
  for c in vecu_cache(integrator.cache)
    resize!(c,i)
  end
  for c in dual_cache(integrator.cache)
    resize!(c.du,i)
    resize!(c.dual_du,i)
  end
  if alg_autodiff(integrator.alg)
    cache.adf = autodiff_setup(cache.rhs,cache.uhold,integrator.alg)
  end
end

function deleteat!(integrator::DDEIntegrator,i::Int)
  for c in full_cache(integrator)
    deleteat!(c,i)
  end
end

function terminate!(integrator::DDEIntegrator)
  integrator.opts.tstops.valtree = typeof(integrator.opts.tstops.valtree)()
end
