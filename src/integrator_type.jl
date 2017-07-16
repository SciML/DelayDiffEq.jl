mutable struct DDEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType,tType,absType,relType,
                             residType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,
                             ProgressType,CacheType,IType,ProbType,NType,O,SType} <:
                                 AbstractDDEIntegrator

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
    fixedpoint_abstol::absType
    fixedpoint_reltol::relType
    resid::residType # This would have to resize for resizing DDE to work
    fixedpoint_norm::NType
    max_fixedpoint_iters::Int
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
    save_idxs::SType
    saveat::Vector{tType}
    save_everystep::Bool
    save_start::Bool
    fsalfirst::rateType
    fsallast::rateType

    # incomplete initialization without fsalfirst and fsallast
    function DDEIntegrator{algType,uType,tType,absType,relType,residType,tTypeNoUnits,
                           tdirType,ksEltype,SolType,rateType,F,ProgressType,CacheType,
                           IType,ProbType,NType,O,SType}(
                               sol,prob,u,k,t,dt,f,uprev,tprev,u_cache,fixedpoint_abstol,
                               fixedpoint_reltol,resid,fixedpoint_norm,max_fixedpoint_iters,
                               alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,
                               dtpropose,tdir,EEst,qold,q11,iter,saveiter,saveiter_dense,
                               prog,cache,kshortsize,just_hit_tstop,accept_step,isout,
                               reeval_fsal,u_modified,opts,integrator,save_idxs,saveat,
                               save_everystep,save_start) where
        {algType<:OrdinaryDiffEqAlgorithm,uType,tType,absType,relType,residType,
         tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,ProgressType,CacheType,IType,
         ProbType,NType,O,SType}

        new(sol,prob,u,k,t,dt,f,uprev,tprev,u_cache,fixedpoint_abstol,fixedpoint_reltol,
            resid,fixedpoint_norm,max_fixedpoint_iters,alg,rate_prototype,notsaveat_idxs,
            dtcache,dtchangeable,dtpropose,tdir,EEst,qold,q11,iter,saveiter,saveiter_dense,
            prog,cache,kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,
            opts,integrator,save_idxs,saveat,save_everystep,save_start)
    end
end

function (integrator::DDEIntegrator)(t, deriv::Type=Val{0}; idxs=nothing)
    OrdinaryDiffEq.current_interpolant(t, integrator, idxs, deriv)
end

function (integrator::DDEIntegrator)(val::AbstractArray, t::Union{Number,AbstractArray},
                                     deriv::Type=Val{0}; idxs=nothing)
    OrdinaryDiffEq.current_interpolant!(val, t, integrator, idxs, deriv)
end
