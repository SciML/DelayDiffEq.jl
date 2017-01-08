type DDEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType<:Union{AbstractArray,Number},tType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,ProgressType,CacheType,IType,ProbType,O} <: AbstractODEIntegrator
  sol::SolType
  prob::ProbType
  u::uType
  k::ksEltype
  t::tType
  dt::tType
  f::F
  uprev::uType
  tprev::tType
  alg::algType
  rate_prototype::rateType
  notsaveat_idxs::Vector{Int}
  dtcache::tType
  dtchangeable::Bool
  dtpropose::tType
  dt_mod::tTypeNoUnits
  tdir::tdirType
  EEst::tTypeNoUnits
  qold::tTypeNoUnits
  q11::tTypeNoUnits
  iter::Int
  saveiter::Int
  saveiter_dense::Int
  prog::ProgressType
  cache::CacheType
  kshortsize::Int
  just_hit_tstop::Bool
  accept_step::Bool
  reeval_fsal::Bool
  u_modified::Bool
  opts::O
  integrator::IType
  fsalfirst::rateType
  fsallast::rateType

  DDEIntegrator(sol,prob,u,k,t,dt,f,uprev,tprev,
      alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,dtpropose,dt_mod,tdir,
      EEst,qold,q11,iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,reeval_fsal,u_modified,integrator,opts) = new(
      sol,prob,u,k,t,dt,f,uprev,tprev,
      alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,dtpropose,dt_mod,tdir,
      EEst,qold,q11,iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,reeval_fsal,u_modified,integrator,opts) # Leave off fsalfirst and last
end

function savevalues!(integrator::DDEIntegrator)
  integrator.integrator.u = integrator.u
  integrator.integrator.k = integrator.k
  integrator.integrator.t = integrator.t
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
  integrator.integrator.fsalfirst = integrator.fsalfirst
  integrator.integrator.t = integrator.t
  integrator.integrator.dt = integrator.dt
  perform_step!(integrator.integrator,integrator.cache,integrator.f)
  integrator.u = integrator.integrator.u
  integrator.fsallast = integrator.integrator.fsallast
  if integrator.opts.adaptive
    integrator.EEst = integrator.integrator.EEst
  end
end

function initialize!(dde_int::DDEIntegrator)
  initialize!(dde_int.integrator,dde_int.cache,dde_int.f)
  dde_int.kshortsize = dde_int.integrator.kshortsize
  dde_int.k = dde_int.integrator.k
  if OrdinaryDiffEq.isfsal(dde_int.alg)
    dde_int.fsalfirst = dde_int.integrator.fsalfirst
  end
end

(integrator::DDEIntegrator)(t) = current_interpolant(t,integrator)
