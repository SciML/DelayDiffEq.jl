function savevalues!(integrator::DDEIntegrator,force_save=false)
  # update ODE integrator with values of DDE integrator
  integrator.integrator.u = integrator.u
  integrator.integrator.k = integrator.k
  integrator.integrator.t = integrator.t

  # add steps for interpolation when needed
  OrdinaryDiffEq.ode_addsteps!(integrator.integrator,integrator.f)

  # update solution of ODE integrator
  savevalues!(integrator.integrator,force_save)
end

function postamble!(integrator::DDEIntegrator)
  # update ODE integrator with values of DDE integrator
  integrator.integrator.u = integrator.u
  integrator.integrator.k = integrator.k
  integrator.integrator.t = integrator.t

  # clean up solution of ODE integrator
  OrdinaryDiffEq.postamble!(integrator.integrator)
end

function perform_step!(integrator::DDEIntegrator)
  integrator.tprev = integrator.t # this is necessary to extrapolate from the current interval

  # update ODE integrator with values of DDE integrator
  integrator.integrator.uprev = integrator.uprev
  integrator.integrator.tprev = integrator.tprev
  integrator.integrator.fsalfirst = integrator.fsalfirst
  integrator.integrator.t = integrator.t
  integrator.integrator.dt = integrator.dt

  # if dt>min lag, then it's explicit so use fixed-point iteration
  if integrator.dt >minimum(integrator.prob.lags)
    # this is done to correct the extrapolation
    t_cache = integrator.t
    tprev_cache = integrator.tprev
    if typeof(integrator.uprev_cache) <: AbstractArray
      copy!(integrator.uprev_cache,integrator.uprev)
    else
      integrator.uprev_cache = integrator.uprev
    end

    if eltype(integrator.u) <: Real # not possible to use nlsolve with units (yet)
      integrator.first_iteration = true # first iteration step needs special updates

      # execute Anderson acceleration of fixed-point iteration
      if typeof(integrator.u) <: Vector
        nlsolve(integrator.iterator,integrator.u;
	        method=:anderson,m=integrator.m,iterations=integrator.max_fixedpoint_iters,xtol=0,ftol=1)
      else
        nlsolve(integrator.iterator,vec(collect(integrator.u));
	        method=:anderson,m=integrator.m,iterations=integrator.max_fixedpoint_iters,xtol=0,ftol=1)
      end
    else # simple Picard iteration for units
      for i in 1:integrator.max_fixedpoint_iters
        # update cache of u
        if typeof(integrator.u) <: AbstractArray
          copy!(integrator.u_cache,integrator.u)
        else
          integrator.u_cache = integrator.u
        end

        perform_step!(integrator,integrator.cache)

	# calculate residuals
        if typeof(integrator.resid) <: AbstractArray
          @. integrator.resid = (integrator.u - integrator.u_cache)/
            @muladd(integrator.fixedpoint_abstol + max(abs(integrator.u),abs(integrator.u_cache))*integrator.fixedpoint_reltol)
        else
          integrator.resid = @. (integrator.u - integrator.u_cache)/
            (integrator.fixedpoint_abstol + max(abs(integrator.u),abs(integrator.u_cache))*integrator.fixedpoint_reltol)
        end

	# special updates of ODE integrator in first iteration step
        if i == 1
          integrator.integrator.tprev = integrator.t
          integrator.integrator.t = integrator.t+integrator.dt
          integrator.integrator.uprev = integrator.u
        end

	# update ODE integrator with values of DDE integrator
        integrator.integrator.u = integrator.u
        integrator.integrator.k = integrator.k

	# stop fixed-point iteration when residuals are small
	integrator.picardnorm(integrator.resid) < 1 && break
      end
    end
    
    # reset to cached values
    integrator.t = t_cache
    integrator.tprev = tprev_cache
    if typeof(integrator.uprev) <: AbstractArray
      copy!(integrator.uprev,integrator.uprev_cache)
    else
      integrator.uprev = integrator.uprev_cache
    end
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

@inline function u_modified!(integrator::DDEIntegrator,bool::Bool)
  integrator.u_modified = bool
end

@inline get_proposed_dt(integrator::DDEIntegrator) = integrator.dtpropose
@inline set_proposed_dt!(integrator::DDEIntegrator,dt) = (integrator.dtpropose = dt)

user_cache(integrator::DDEIntegrator) = user_cache(integrator)
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
