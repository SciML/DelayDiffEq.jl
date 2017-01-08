function init{uType,tType,isinplace,algType<:AbstractMethodOfStepsAlgorithm,lType,F,H}(
  prob::AbstractDDEProblem{uType,tType,lType,isinplace,F,H},
  alg::algType,timeseries_init=uType[],ts_init=tType[],ks_init=[];
  d_discontinuities = tType[],
  dtmax=tType((prob.tspan[end]-prob.tspan[1])),
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

  if dt == zero(dt) && adaptive
    dt = tType(ode_determine_initdt(prob.u0,prob.tspan[1],integrator.tdir,dtmax,integrator.abstol,integrator.reltol,integrator.internalnorm,dde_f,OrdinaryDiffEq.alg_order(order)))
  end
  integrator.dt = dt





  dde_int = DDEIntegrator{typeof(integrator.alg),
                             uType,tType,
                             tTypeNoUnits,typeof(integrator.tdir),
                             typeof(integrator.k),typeof(integrator.sol),
                             typeof(integrator.rate_prototype),
                             typeof(dde_f),typeof(integrator.prog),
                             typeof(integrator.cache),
                             typeof(integrator),typeof(prob),
                             typeof(integrator.opts)}(
      integrator.sol,prob,integrator.u,integrator.k,integrator.t,integrator.dt,
      dde_f,integrator.uprev,integrator.tprev,
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
