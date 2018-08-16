mutable struct DDEIntegrator{algType<:OrdinaryDiffEq.OrdinaryDiffEqAlgorithm,
                             uType,tType,P,eigenType,absType,relType,
                             residType,tTypeNoUnits,tdirType,ksEltype,
                             SolType,F,CacheType,
                             IType,NType,O,tstopsType,
                             FSALType} <: AbstractDDEIntegrator
    sol::SolType
    u::uType
    k::ksEltype
    t::tType
    dt::tType
    f::F
    p::P
    uprev::uType
    uprev2::uType
    tprev::tType
    prev_idx::Int
    prev2_idx::Int
    fixedpoint_abstol::absType
    fixedpoint_reltol::relType
    resid::residType # This would have to resize for resizing DDE to work
    fixedpoint_norm::NType
    max_fixedpoint_iters::Int
    saveat::tstopsType
    tracked_discontinuities::Vector{Discontinuity{tType}}
    alg::algType
    dtcache::tType
    dtchangeable::Bool
    dtpropose::tType
    tdir::tdirType
    eigen_est::eigenType
    EEst::tTypeNoUnits
    qold::tTypeNoUnits
    q11::tTypeNoUnits
    erracc::tTypeNoUnits
    dtacc::tType
    success_iter::Int
    iter::Int
    saveiter::Int
    saveiter_dense::Int
    cache::CacheType
    kshortsize::Int
    force_stepfail::Bool
    last_stepfail::Bool
    just_hit_tstop::Bool
    event_last_time::Int
    accept_step::Bool
    isout::Bool
    reeval_fsal::Bool
    u_modified::Bool
    opts::O
    integrator::IType
    fsalfirst::FSALType
    fsallast::FSALType

    # incomplete initialization without fsalfirst and fsallast
    function DDEIntegrator{algType,uType,tType,P,eigenType,absType,relType,
                           residType,tTypeNoUnits,
                           tdirType,ksEltype,SolType,F,CacheType,IType,
                           NType,O,tstopsType,FSALType}(
                               sol,u,k,t,dt,f,p,uprev,uprev2,tprev,prev_idx,prev2_idx,
                               fixedpoint_abstol,fixedpoint_reltol,resid,fixedpoint_norm,
                               max_fixedpoint_iters,saveat,tracked_discontinuities,alg,
                               dtcache,dtchangeable,dtpropose,tdir,eigen_est,EEst,qold,
                               q11,erracc,dtacc,success_iter,iter,saveiter,saveiter_dense,
                               cache,kshortsize,force_stepfail,last_stepfail,
                               just_hit_tstop,event_last_time,accept_step,isout,
                               reeval_fsal,u_modified,opts,
                               integrator) where
        {algType,uType,tType,P,eigenType,absType,relType,residType,tTypeNoUnits,
        tdirType,ksEltype,
         SolType,F,CacheType,IType,NType,O,tstopsType,FSALType}

        new(sol,u,k,t,dt,f,p,uprev,uprev2,tprev,prev_idx,prev2_idx,fixedpoint_abstol,
            fixedpoint_reltol,resid,fixedpoint_norm,max_fixedpoint_iters,saveat,
            tracked_discontinuities,alg,dtcache,dtchangeable,dtpropose,tdir,eigen_est,
            EEst,qold,q11,erracc,dtacc,success_iter,iter,saveiter,saveiter_dense,cache,
            kshortsize,force_stepfail,last_stepfail,just_hit_tstop,accept_step,
            event_last_time,isout,reeval_fsal,u_modified,opts,integrator)
    end
end

function (integrator::DDEIntegrator)(t, deriv::Type=Val{0}; idxs=nothing)
    OrdinaryDiffEq.current_interpolant(t, integrator, idxs, deriv)
end

function (integrator::DDEIntegrator)(val::AbstractArray, t::Union{Number,AbstractArray},
                                     deriv::Type=Val{0}; idxs=nothing)
    OrdinaryDiffEq.current_interpolant!(val, t, integrator, idxs, deriv)
end
