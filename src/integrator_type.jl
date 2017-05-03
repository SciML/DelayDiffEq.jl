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
  picardabstol::absType
  picardreltol::relType
  resid::residType # This would have to resize for resizing DDE to work
  picardnorm::NType
  max_picard_iters::Int
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

  (::Type{DDEIntegrator{algType,uType,tType,absType,relType,
                   residType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,
                   ProgressType,CacheType,IType,ProbType,NType,O}}){
                   algType,uType,tType,absType,
                   relType,residType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,ProgressType,
                  CacheType,IType,ProbType,NType,O}(sol,prob,u,k,t,dt,f,uprev,tprev,u_cache,
      picardabstol,picardreltol,resid,picardnorm,max_picard_iters,
      alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,dtpropose,tdir,
      EEst,qold,q11,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,
      integrator,opts) = new{algType,uType,tType,absType,relType,residType,tTypeNoUnits,tdirType,
      ksEltype,SolType,rateType,F,ProgressType,CacheType,IType,ProbType,NType,O}(
      sol,prob,u,k,t,dt,f,uprev,tprev,u_cache,
      picardabstol,picardreltol,resid,picardnorm,max_picard_iters,
      alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,dtpropose,tdir,
      EEst,qold,q11,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,integrator,opts) # Leave off fsalfirst and last
end

function (integrator::DDEIntegrator)(t,deriv::Type=Val{0};idxs=size(integrator.uprev))
  OrdinaryDiffEq.current_interpolant(t,integrator,idxs,deriv)
end
(integrator::DDEIntegrator)(val::AbstractArray,t::Union{Number,AbstractArray},deriv::Type=Val{0};idxs=eachindex(integrator.uprev)) = OrdinaryDiffEq.current_interpolant!(val,t,integrator,idxs,deriv)
