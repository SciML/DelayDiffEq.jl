type DDEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType,tType,absType,relType,residType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,ProgressType,CacheType,IType,ProbType,NType,O} <: AbstractDDEIntegrator
  sol::SolType
  prob::ProbType
  u::uType
  k::ksEltype
  t::tType
  dt::tType
  f::F
  uprev::uType
  tprev::tType
  u_cache::uType
  uprev_cache::uType
  fixedpoint_abstol::absType
  fixedpoint_reltol::relType
  resid::residType # This would have to resize for resizing DDE to work
  picardnorm::NType
  max_fixedpoint_iters::Int
  m::Int
  alg::algType
  rate_prototype::rateType
  notsaveat_idxs::Vector{Int}
  dtcache::tType
  dtchangeable::Bool
  dtpropose::tType
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
  isout::Bool
  reeval_fsal::Bool
  u_modified::Bool
  opts::O
  integrator::IType
  fsalfirst::rateType
  fsallast::rateType
  first_iteration::Bool
  iterator::IterationFunction{DDEIntegrator{algType,uType,tType,absType,relType,residType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,ProgressType,CacheType,IType,ProbType,NType,O}}

  (::Type{DDEIntegrator{algType,uType,tType,absType,relType,
                   residType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,
                   ProgressType,CacheType,IType,ProbType,NType,O}}){
                   algType,uType,tType,absType,
                   relType,residType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,ProgressType,
                  CacheType,IType,ProbType,NType,O}(sol,prob,u,k,t,dt,f,uprev,tprev,u_cache,uprev_cache,
      fixedpoint_abstol,fixedpoint_reltol,resid,picardnorm,max_fixedpoint_iters,m,
      alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,dtpropose,tdir,
      EEst,qold,q11,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,
      opts,integrator) = begin
    dde_int = new{algType,uType,tType,absType,relType,residType,tTypeNoUnits,tdirType,
      ksEltype,SolType,rateType,F,ProgressType,CacheType,IType,ProbType,NType,O}(
      sol,prob,u,k,t,dt,f,uprev,tprev,u_cache,uprev_cache,
      fixedpoint_abstol,fixedpoint_reltol,resid,picardnorm,max_fixedpoint_iters,m,
      alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,dtpropose,tdir,
      EEst,qold,q11,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,opts,integrator) # Leave off fsalfirst, fsallast, first_iteration and iterator
    dde_int.iterator = IterationFunction{DDEIntegrator{algType,uType,tType,absType,relType,residType,tTypeNoUnits,tdirType,
                                                      ksEltype,SolType,rateType,F,ProgressType,CacheType,IType,ProbType,NType,O}}(dde_int)
    dde_int
  end
end

function (integrator::DDEIntegrator)(t,deriv::Type=Val{0};idxs=nothing)
  OrdinaryDiffEq.current_interpolant(t,integrator,idxs,deriv)
end
(integrator::DDEIntegrator)(val::AbstractArray,t::Union{Number,AbstractArray},deriv::Type=Val{0};idxs=nothing) = OrdinaryDiffEq.current_interpolant!(val,t,integrator,idxs,deriv)
