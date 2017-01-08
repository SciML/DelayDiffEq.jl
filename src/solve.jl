function init{uType,tType,isinplace,algType<:AbstractMethodOfStepsAlgorithm,lType,F,H}(
  prob::AbstractDDEProblem{uType,tType,lType,isinplace,F,H},
  alg::algType,timeseries_init=uType[],ts_init=tType[],ks_init=[];
  d_discontinuities = tType[],
  dtmax=tType(10*minimum(prob.lags)),
  dt = tType(0),
  kwargs...)

  # Add to the discontinuties vector the lag locations
  d_discontinuities = [d_discontinuities;vec([prob.tspan[1]+i*τ for i=1:alg_order(alg),τ in prob.lags])]

  # If it's constrained, then no Picard iteration, and thus `dtmax` should match max lag size
  if isconstrained(alg)
    dtmax = min(dtmax,prob.lags...)
  end

  tTypeNoUnits   = typeof(recursive_one(prob.tspan[1]))

  # Bootstrap the Integrator Using An ODEProblem
  if typeof(prob) <: AbstractDDEProblem
    ode_prob = ODEProblem(prob.f,prob.u0,prob.tspan;iip=isinplace)
  elseif typeof(prob) <: AbstractDDETestProblem
    ode_prob = ODETestProblem(prob.f,prob.u0,prob.analytic,prob.tspan;iip=isinplace)
  end
  integrator = init(ode_prob,alg.alg;dt=1,initialize_integrator=false,
                    d_discontinuities=d_discontinuities,
                    dtmax=dtmax,
                    kwargs...)
  h = HistoryFunction(prob.h,integrator.sol,integrator)
  if isinplace
    dde_f = (t,u,du) -> prob.f(t,u,h,du)
  else
    dde_f = (t,u) -> prob.f(t,u,h)
  end

  if typeof(alg.alg) <: OrdinaryDiffEqCompositeAlgorithm
    id = OrdinaryDiffEq.CompositeInterpolationData(integrator.sol.interp,dde_f)
  else
    id = OrdinaryDiffEq.InterpolationData(integrator.sol.interp,dde_f)
  end

  if typeof(alg.alg) <: OrdinaryDiffEqCompositeAlgorithm
    sol = build_solution(integrator.sol.prob,
                         integrator.sol.alg,
                         integrator.sol.t,
                         integrator.sol.u,
                         dense=integrator.sol.dense,
                         k=integrator.sol.k,
                         interp=id,
                         alg_choice=integrator.sol.alg_choice,
                         calculate_error = false)
  else
    sol = build_solution(integrator.sol.prob,
                         integrator.sol.alg,
                         integrator.sol.t,
                         integrator.sol.u,
                         dense=integrator.sol.dense,
                         k=integrator.sol.k,
                         interp=id,
                         calculate_error = false)
  end

  h2 = HistoryFunction(prob.h,sol,integrator)
  if isinplace
    dde_f2 = (t,u,du) -> prob.f(t,u,h2,du)
  else
    dde_f2 = (t,u) -> prob.f(t,u,h2)
  end

  if dt == zero(dt) && integrator.opts.adaptive
    dt = tType(OrdinaryDiffEq.ode_determine_initdt(prob.u0,prob.tspan[1],
              integrator.tdir,dtmax,integrator.opts.abstol,
              integrator.opts.reltol,integrator.opts.internalnorm,
              dde_f,OrdinaryDiffEq.alg_order(alg)))
  end
  integrator.dt = dt

  if typeof(alg.picardabstol) <: Void
    picardabstol_internal = integrator.opts.abstol
  else
    picardabstol_internal = alg.picardabstol
  end
  if typeof(alg.picardreltol) <: Void
    picardreltol_internal = integrator.opts.reltol
  else
    picardreltol_internal = alg.picardreltol
  end
  if typeof(alg.picardnorm) <: Void
    picardnorm = integrator.opts.internalnorm
  end


  uEltypeNoUnits = typeof(recursive_one(integrator.u))

  if typeof(integrator.u) <: AbstractArray
    resid = similar(integrator.u,uEltypeNoUnits)
    u_cache = similar(integrator.u)
  else
    resid = uEltypeNoUnits(1)
    u_cache = one(uType)
  end



  dde_int = DDEIntegrator{typeof(integrator.alg),
                             uType,tType,
                             typeof(picardabstol_internal),
                             typeof(picardreltol_internal),
                             typeof(resid),
                             tTypeNoUnits,typeof(integrator.tdir),
                             typeof(integrator.k),typeof(sol),
                             typeof(integrator.rate_prototype),
                             typeof(dde_f2),typeof(integrator.prog),
                             typeof(integrator.cache),
                             typeof(integrator),typeof(prob),
                             typeof(picardnorm),
                             typeof(integrator.opts)}(
      sol,prob,integrator.u,integrator.k,integrator.t,integrator.dt,
      dde_f2,integrator.uprev,integrator.tprev,u_cache,
      eltype(integrator.u)(picardabstol_internal),uEltypeNoUnits(picardreltol_internal),
      resid,picardnorm,alg.max_picard_iters,
      integrator.alg,integrator.rate_prototype,integrator.notsaveat_idxs,integrator.dtcache,
      integrator.dtchangeable,integrator.dtpropose,integrator.dt_mod,integrator.tdir,
      integrator.EEst,integrator.qold,integrator.q11,integrator.iter,integrator.saveiter,
      integrator.saveiter_dense,integrator.prog,integrator.cache,
      integrator.kshortsize,integrator.just_hit_tstop,integrator.accept_step,
      integrator.reeval_fsal,integrator.u_modified,integrator.opts,integrator) # Leave off fsalfirst and last

  initialize!(dde_int)
  dde_int
end

function solve!(dde_int::DDEIntegrator)
  @inbounds while !isempty(dde_int.opts.tstops)
    while dde_int.tdir*dde_int.t < dde_int.tdir*top(dde_int.opts.tstops)
      loopheader!(dde_int)
      perform_step!(dde_int)
      loopfooter!(dde_int)
      if isempty(dde_int.opts.tstops)
        break
      end
    end
    handle_tstop!(dde_int)
  end

  postamble!(dde_int)
  if typeof(dde_int.prob) <: AbstractDDETestProblem
    u_analytic = [dde_int.prob.analytic(t,dde_int.sol[1]) for t in dde_int.sol.t]
    errors = Dict{Symbol,eltype(dde_int.u)}()
    sol = build_solution(dde_int.sol::AbstractODESolution,u_analytic,errors)
    calculate_solution_errors!(sol;fill_uanalytic=false,timeseries_errors=dde_int.opts.timeseries_errors,dense_errors=dde_int.opts.dense_errors)
    return sol
  else
    return dde_int.sol
  end
end

function solve{uType,tType,isinplace,algType<:AbstractMethodOfStepsAlgorithm,lType,F,H}(
  prob::AbstractDDEProblem{uType,tType,lType,isinplace,F,H},
  alg::algType,timeseries_init=uType[],ts_init=tType[],ks_init=[];kwargs...)

  integrator = init(prob,alg,timeseries_init,ts_init,ks_init;kwargs...)
  solve!(integrator)
end
