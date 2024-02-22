mutable struct HistoryODEIntegrator{
    algType, IIP, uType, tType, tdirType, ksEltype, SolType,
    CacheType, DV} <:
               AbstractODEIntegrator{algType, IIP, uType, tType}
    sol::SolType
    u::uType
    k::ksEltype
    t::tType
    dt::tType
    uprev::uType
    tprev::tType
    alg::algType
    dtcache::tType
    tdir::tdirType
    saveiter::Int
    saveiter_dense::Int
    cache::CacheType
    differential_vars::DV
end

function (integrator::HistoryODEIntegrator)(t, deriv::Type = Val{0}; idxs = nothing)
    OrdinaryDiffEq.current_interpolant(t, integrator, idxs, deriv)
end

function (integrator::HistoryODEIntegrator)(val::AbstractArray,
        t::Union{Number, AbstractArray},
        deriv::Type = Val{0}; idxs = nothing)
    OrdinaryDiffEq.current_interpolant!(val, t, integrator, idxs, deriv)
end

mutable struct DDEIntegrator{algType, IIP, uType, tType, P, eigenType, tTypeNoUnits,
    tdirType,
    ksEltype, SolType, F, CacheType, IType, FP, O, dAbsType,
    dRelType, H,
    tstopsType, discType, FSALType, EventErrorType,
    CallbackCacheType, DV} <:
               AbstractDDEIntegrator{algType, IIP, uType, tType}
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
    fpsolver::FP
    order_discontinuity_t0::Int
    "Discontinuities tracked by callback."
    tracked_discontinuities::Vector{Discontinuity{tType, Int}}
    discontinuity_interp_points::Int
    discontinuity_abstol::dAbsType
    discontinuity_reltol::dRelType
    "Future time stops for propagated discontinuities."
    tstops_propagated::tstopsType
    "Future propagated discontinuities."
    d_discontinuities_propagated::discType
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
    callback_cache::CallbackCacheType
    kshortsize::Int
    force_stepfail::Bool
    last_stepfail::Bool
    just_hit_tstop::Bool
    do_error_check::Bool
    event_last_time::Int
    vector_event_last_time::Int
    last_event_error::EventErrorType
    accept_step::Bool
    isout::Bool
    reeval_fsal::Bool
    u_modified::Bool
    isdae::Bool
    opts::O
    stats::SciMLBase.DEStats
    history::H
    differential_vars::DV
    integrator::IType
    fsalfirst::FSALType
    fsallast::FSALType

    # incomplete initialization without fsalfirst and fsallast
    function DDEIntegrator{algType, IIP, uType, tType, P, eigenType, tTypeNoUnits,
            tdirType, ksEltype, SolType, F, CacheType, IType, FP,
            O, dAbsType, dRelType, H, tstopsType, discType,
            FSALType, EventErrorType, CallbackCacheType, DV}(sol, u, k, t, dt, f,
            p, uprev, uprev2,
            tprev, prev_idx,
            prev2_idx,
            fpsolver,
            order_discontinuity_t0,
            tracked_discontinuities,
            discontinuity_interp_points,
            discontinuity_abstol,
            discontinuity_reltol,
            tstops_propagated,
            d_discontinuities_propagated,
            alg, dtcache,
            dtchangeable,
            dtpropose, tdir,
            eigen_est, EEst,
            qold,
            q11, erracc, dtacc,
            success_iter, iter,
            saveiter,
            saveiter_dense,
            cache,
            callback_cache,
            kshortsize,
            force_stepfail,
            last_stepfail,
            just_hit_tstop,
            do_error_check,
            event_last_time,
            vector_event_last_time,
            last_event_error,
            accept_step, isout,
            reeval_fsal,
            u_modified, isdae,
            opts, stats,
            history,
            differential_vars,
            integrator) where
            {algType, IIP, uType, tType, P, eigenType, tTypeNoUnits, tdirType, ksEltype,
            SolType, F,
            CacheType, IType, FP, O, dAbsType, dRelType, H, tstopsType, discType,
            FSALType, EventErrorType, CallbackCacheType, DV}
        new{algType, IIP, uType, tType, P, eigenType, tTypeNoUnits, tdirType, ksEltype,
            SolType, F,
            CacheType, IType, FP, O, dAbsType, dRelType, H, tstopsType, discType, FSALType,
            EventErrorType, CallbackCacheType, DV}(
            sol, u, k, t, dt, f, p, uprev, uprev2, tprev,
            prev_idx, prev2_idx, fpsolver,
            order_discontinuity_t0,
            tracked_discontinuities,
            discontinuity_interp_points,
            discontinuity_abstol, discontinuity_reltol,
            tstops_propagated,
            d_discontinuities_propagated, alg, dtcache,
            dtchangeable, dtpropose, tdir,
            eigen_est, EEst, qold, q11, erracc, dtacc,
            success_iter, iter, saveiter, saveiter_dense,
            cache, callback_cache, kshortsize,
            force_stepfail, last_stepfail,
            just_hit_tstop,
            do_error_check, event_last_time,
            vector_event_last_time,
            last_event_error, accept_step, isout,
            reeval_fsal, u_modified, isdae, opts,
            stats, history, differential_vars, integrator)
    end
end

function (integrator::DDEIntegrator)(t, deriv::Type = Val{0}; idxs = nothing)
    OrdinaryDiffEq.current_interpolant(t, integrator, idxs, deriv)
end

function (integrator::DDEIntegrator)(val::AbstractArray, t::Union{Number, AbstractArray},
        deriv::Type = Val{0}; idxs = nothing)
    OrdinaryDiffEq.current_interpolant!(val, t, integrator, idxs, deriv)
end
