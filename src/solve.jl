function init(prob::AbstractDDEProblem{uType,tType,lType,isinplace}, alg::algType,
              timeseries_init=uType[], ts_init=tType[], ks_init=[];
              d_discontinuities=tType[], dtmax=tType(7*minimum(prob.lags)), dt=zero(tType),
              saveat=tType[], tstops=tType[], save_idxs=nothing,
              save_everystep=isempty(saveat), save_start=true, kwargs...) where
    {uType,tType,isinplace,algType<:AbstractMethodOfStepsAlgorithm,lType}

    # add lag locations to discontinuities vector
    d_discontinuities = [d_discontinuities; compute_discontinuity_tree(prob.lags, alg,
                                                                       prob.tspan[1])]

    # no fixed-point iterations for constrained algorithms,
    # and thus `dtmax` should match minimal lag
    if isconstrained(alg)
        dtmax = min(dtmax,prob.lags...)
    end

    # bootstrap the integrator using an ODE problem, but do not initialize it since
    # ODE solvers only accept functions f(t,u,du) or f(t,u) without history function
    # save every step and every index in the solution of the integrator, since they are
    # required for evaluation of past values during integration
    ode_prob = ODEProblem(prob.f, prob.u0, prob.tspan; iip=isinplace)
    integrator = init(ode_prob, alg.alg; dt=1, initialize_integrator=false,
                      d_discontinuities=d_discontinuities, dtmax=dtmax, kwargs...)

    # create new solution based on this integrator with an interpolation function of the
    # expected form f(t,u,du) or f(t,u) which already includes information about the
    # history function of the DDE problem, the current solution of the integrator, and
    # the extrapolation of the integrator for the future
    interp_h = HistoryFunction(prob.h, integrator.sol, integrator)
    if isinplace
        interp_f = (t,u,du) -> prob.f(t,u,interp_h,du)
    else
        interp_f = (t,u) -> prob.f(t,u,interp_h)
    end

    if typeof(alg.alg) <: OrdinaryDiffEqCompositeAlgorithm
        ode_interp_data = OrdinaryDiffEq.CompositeInterpolationData(integrator.sol.interp,
                                                                    interp_f)
    else
        ode_interp_data = OrdinaryDiffEq.InterpolationData(integrator.sol.interp, interp_f)
    end

    if typeof(alg.alg) <: OrdinaryDiffEqCompositeAlgorithm
        ode_sol = build_solution(prob, integrator.sol.alg, integrator.sol.t,
                                 integrator.sol.u, dense=integrator.sol.dense,
                                 k=integrator.sol.k, interp=ode_interp_data,
                                 alg_choice=integrator.sol.alg_choice, calculate_error=false)
    else
        ode_sol = build_solution(prob, integrator.sol.alg, integrator.sol.t,
                                 integrator.sol.u, dense=integrator.sol.dense,
                                 k=integrator.sol.k, interp=ode_interp_data,
                                 calculate_error=false)
    end

    # use this improved solution together with the given history function and the integrator
    # to create a problem function of the DDE with all available history information that is
    # of the form f(t,u,du) or f(t,u) such that ODE algorithms can be applied
    dde_h = HistoryFunction(prob.h, ode_sol, integrator)
    if isinplace
        dde_f = (t,u,du) -> prob.f(t,u,dde_h,du)
    else
        dde_f = (t,u) -> prob.f(t,u,dde_h)
    end

    # use ODE problem to initialize prototype of DDE integrator
    ode_prob = ODEProblem(dde_f, prob.u0, prob.tspan; iip=isinplace)
    proto_int = init(ode_prob, alg.alg; dt=1, initialize_integrator=false,
                     d_discontinuities=d_discontinuities, saveat=saveat, dtmax=dtmax,
                     save_idxs=save_idxs, save_everystep=save_everystep,
                     save_start=save_start, kwargs...)

    # if time step is not set and prototype uses adaptive step sizes,
    # let function with history support calculate initial time step
    if dt == zero(dt) && proto_int.opts.adaptive
        dt = tType(OrdinaryDiffEq.ode_determine_initdt(prob.u0, prob.tspan[1],
                                                       proto_int.tdir, minimum(prob.lags),
                                                       proto_int.opts.abstol,
                                                       proto_int.opts.reltol,
                                                       proto_int.opts.internalnorm,
                                                       ode_prob,
                                                       OrdinaryDiffEq.alg_order(alg)))
    end
    proto_int.dt = dt

    # create solution for DDE integrator, using solution of prototype and cache of ODE
    # integrator
    if typeof(alg.alg) <: OrdinaryDiffEqCompositeAlgorithm
        CompositeInterpolationData(f,timeseries,ts,ks,alg_choice,notsaveat_idxs,dense,cache)
        dde_interp_data = OrdinaryDiffEq.CompositeInterpolationData(
            dde_f, proto_int.sol.interp.timeseries, proto_int.sol.interp.ts,
            proto_int.sol.interp.ks, proto_int.sol.interp.alg_choice,
            proto_int.sol.interp.notsaveat_idxs, proto_int.sol.interp.dense,
            integrator.cache)
    else
        dde_interp_data = OrdinaryDiffEq.InterpolationData(
            dde_f, proto_int.sol.interp.timeseries, proto_int.sol.interp.ts,
            proto_int.sol.interp.ks, proto_int.sol.interp.notsaveat_idxs,
            proto_int.sol.interp.dense, integrator.cache)
    end

    if typeof(alg.alg) <: OrdinaryDiffEqCompositeAlgorithm
        dde_sol = build_solution(prob, proto_int.sol.alg, proto_int.sol.t, proto_int.sol.u,
                                 dense=proto_int.sol.dense, k=proto_int.sol.k,
                                 interp=dde_interp_data,
                                 alg_choice=proto_int.alg_choice, calculate_error=false)
    else
        dde_sol = build_solution(prob, proto_int.sol.alg, proto_int.sol.t, proto_int.sol.u,
                                 dense=proto_int.sol.dense, k=proto_int.sol.k,
                                 interp=dde_interp_data, calculate_error=false)
    end

    # absolut tolerance for fixed-point iterations has to be of same type as elements of u
    # in particular important for calculations with units
    if typeof(alg.fixedpoint_abstol) <: Void
        fixedpoint_abstol_internal = map(eltype(uType), proto_int.opts.abstol)
    else
        fixedpoint_abstol_internal = map(eltype(uType), alg.fixedpoint_abstol)
    end

    # use norm of ODE integrator if no norm for fixed-point iterations is specified
    if typeof(alg.fixedpoint_norm) <: Void
        fixedpoint_norm = proto_int.opts.internalnorm
    end

    # derive unitless types
    uEltypeNoUnits = typeof(recursive_one(integrator.u))
    tTypeNoUnits = typeof(recursive_one(prob.tspan[1]))

    # relative tolerance for fixed-point iterations has to be of same type as elements of u
    # without units
    # in particular important for calculations with units
    if typeof(alg.fixedpoint_reltol) <: Void
        fixedpoint_reltol_internal = map(uEltypeNoUnits, proto_int.opts.reltol)
    else
        fixedpoint_reltol_internal = map(uEltypeNoUnits, alg.fixedpoint_reltol)
    end

    # create containers for residuals and to cache u with correct dimensions and types
    # in particular for calculations with units residuals have to be unitless
    if typeof(integrator.u) <: AbstractArray
        resid = similar(integrator.u, uEltypeNoUnits)
        u_cache = similar(integrator.u)
    else
        resid = one(uEltypeNoUnits)
        u_cache = oneunit(eltype(uType))
    end

    # create DDE integrator combining the new defined problem function with history
    # information, the improved solution, parameters of the prototype, some initial values
    # of the ODE integrator, and parameters for fixed-point iterations
    # do not initialize fsalfirst and fsallast
    dde_int = DDEIntegrator{typeof(proto_int.alg),uType,tType,
                            typeof(fixedpoint_abstol_internal),
                            typeof(fixedpoint_reltol_internal),typeof(resid),tTypeNoUnits,
                            typeof(proto_int.tdir),typeof(integrator.k),typeof(dde_sol),
                            typeof(proto_int.rate_prototype),typeof(dde_f),
                            typeof(proto_int.prog),typeof(integrator.cache),
                            typeof(integrator),typeof(prob),typeof(fixedpoint_norm),
                            typeof(proto_int.opts)}(
                                dde_sol, prob, integrator.u, integrator.k, proto_int.t,
                                proto_int.dt, dde_f, integrator.uprev, proto_int.tprev,
                                u_cache, fixedpoint_abstol_internal,
                                fixedpoint_reltol_internal, resid, fixedpoint_norm,
                                alg.max_fixedpoint_iters, proto_int.alg,
                                proto_int.rate_prototype, proto_int.notsaveat_idxs,
                                proto_int.dtcache, proto_int.dtchangeable,
                                proto_int.dtpropose, proto_int.tdir, proto_int.EEst,
                                proto_int.qold, proto_int.q11, proto_int.iter,
                                proto_int.saveiter, proto_int.saveiter_dense,
                                proto_int.prog, integrator.cache, proto_int.kshortsize,
                                proto_int.just_hit_tstop, proto_int.accept_step,
                                proto_int.isout, proto_int.reeval_fsal,
                                proto_int.u_modified, proto_int.opts, integrator)

    # set up additional initial values of newly created DDE integrator
    # (such as fsalfirst) and its callbacks
    initialize!(dde_int)
    initialize!(integrator.opts.callback, proto_int.t, integrator.u, dde_int)

    dde_int
end

function solve!(integrator::DDEIntegrator)
    # step over all stopping time points, similar to solving with ODE integrators
    @inbounds while !isempty(integrator.opts.tstops)
        while integrator.tdir * integrator.t < integrator.tdir * top(integrator.opts.tstops)
            # apply step or adapt step size
            loopheader!(integrator)

            # calculate next step
            perform_step!(integrator)

            # calculate proposed next step size, handle callbacks, and update solution
            loopfooter!(integrator)

            if isempty(integrator.opts.tstops)
                break
            end
        end

        # remove hit or passed stopping time points
        handle_tstop!(integrator)
    end

    # clean up solution
    postamble!(integrator)

    # calculate and add errors to solution if analytic solution to problem exists
    # can not use same mechanism as for ODE problems since analytic solutions for DDE depend
    # on initial value u0
    if has_analytic(integrator.prob.f)
        u_analytic = [integrator.prob.f(Val{:analytic}, t, integrator.sol[1])
                      for t in integrator.sol.t]
        errors = Dict{Symbol,eltype(integrator.u)}()
        sol = build_solution(integrator.sol::AbstractODESolution, u_analytic, errors)
        calculate_solution_errors!(sol; fill_uanalytic=false,
                                   timeseries_errors=integrator.opts.timeseries_errors,
                                   dense_errors=integrator.opts.dense_errors)
        sol.retcode = :Success
        return sol
    else
        integrator.sol.retcode = :Success
        return integrator.sol
    end
end

function solve(prob::AbstractDDEProblem{uType,tType,lType,isinplace}, alg::algType,
               timeseries_init=uType[], ts_init=tType[], ks_init=[]; kwargs...) where
    {uType,tType,isinplace,algType<:AbstractMethodOfStepsAlgorithm,lType}

    integrator = init(prob, alg, timeseries_init, ts_init, ks_init; kwargs...)
    solve!(integrator)
end
