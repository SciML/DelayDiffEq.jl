function init(prob::AbstractDDEProblem{uType,tType,lType,isinplace}, alg::algType,
              timeseries_init=uType[], ts_init=tType[], ks_init=[];
              d_discontinuities=tType[], dtmax=tType(7*minimum(prob.lags)), dt=zero(tType),
              saveat=tType[], save_idxs=nothing, save_everystep=isempty(saveat),
              save_start=true, kwargs...) where
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
        interp_data = OrdinaryDiffEq.CompositeInterpolationData(integrator.sol.interp,
                                                                interp_f)
    else
        interp_data = OrdinaryDiffEq.InterpolationData(integrator.sol.interp,
                                                       interp_f)
    end

    if typeof(alg.alg) <: OrdinaryDiffEqCompositeAlgorithm
        sol = build_solution(prob, integrator.sol.alg, integrator.sol.t, integrator.sol.u,
                             dense=integrator.sol.dense, k=integrator.sol.k,
                             interp=interp_data, alg_choice=integrator.sol.alg_choice,
                             calculate_error = false)
    else
        sol = build_solution(prob, integrator.sol.alg, integrator.sol.t, integrator.sol.u,
                             dense=integrator.sol.dense, k=integrator.sol.k,
                             interp=interp_data, calculate_error = false)
    end

    # use this improved solution together with the given history function and the integrator
    # to create a problem function of the DDE with all available history information that is
    # of the form f(t,u,du) or f(t,u) such that ODE algorithms can be applied
    dde_h = HistoryFunction(prob.h, sol, integrator)
    if isinplace
        dde_f = (t,u,du) -> prob.f(t,u,dde_h,du)
    else
        dde_f = (t,u) -> prob.f(t,u,dde_h)
    end

    # if time step is not set and ODE integrator has adaptive step size,
    # let ODE problem with same parameters and newly created function with history support
    # calculate initial time step
    if dt == zero(dt) && integrator.opts.adaptive
        ode_prob = ODEProblem(dde_f, prob.u0, prob.tspan)
        dt = tType(OrdinaryDiffEq.ode_determine_initdt(prob.u0, prob.tspan[1],
                                                       integrator.tdir, minimum(prob.lags),
                                                       integrator.opts.abstol,
                                                       integrator.opts.reltol,
                                                       integrator.opts.internalnorm,
                                                       ode_prob,
                                                       OrdinaryDiffEq.alg_order(alg)))
    end
    integrator.dt = dt

    # absolut tolerance for fixed-point iterations has to be of same type as elements of u
    # in particular important for calculations with units
    if typeof(alg.fixedpoint_abstol) <: Void
        fixedpoint_abstol_internal = map(eltype(uType), integrator.opts.abstol)
    else
        fixedpoint_abstol_internal = map(eltype(uType), alg.fixedpoint_abstol)
    end

    # use norm of ODE integrator if no norm for fixed-point iterations is specified
    if typeof(alg.fixedpoint_norm) <: Void
        fixedpoint_norm = integrator.opts.internalnorm
    end

    # derive unitless types
    uEltypeNoUnits = typeof(recursive_one(integrator.u))
    tTypeNoUnits = typeof(recursive_one(prob.tspan[1]))

    # relative tolerance for fixed-point iterations has to be of same type as elements of u
    # without units
    # in particular important for calculations with units
    if typeof(alg.fixedpoint_reltol) <: Void
        fixedpoint_reltol_internal = map(uEltypeNoUnits, integrator.opts.reltol)
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

    # create heap of additional time points that will be contained in the solution
    # exclude the end point because of floating point issues and the starting point since it
    # is controlled by save_start
    if typeof(saveat) <: Number
        saveat_internal = collect(tType,
                                  (prob.tspan[1] + integrator.tdir * abs(saveat)):
                                  integrator.tdir * abs(saveat):
                                  (prob.tspan[end] - integrator.tdir * abs(saveat)))
    else
        saveat_internal = collect(tType, Iterators.filter(
            x -> integrator.tdir * prob.tspan[1] < integrator.tdir * x < integrator.tdir *
            prob.tspan[end],
            saveat))
    end

    if integrator.tdir > 0
        heapify!(saveat_internal)
    else
        heapify!(saveat_internal, Base.Order.Reverse)
    end

    # create DDE integrator combining the new defined problem function with history
    # information, the improved solution, the parameters of the ODE integrator, and
    # parameters of fixed-point iteration
    # do not initialize fsalfirst and fsallast
    dde_int = DDEIntegrator{typeof(integrator.alg),uType,tType,
                            typeof(fixedpoint_abstol_internal),
                            typeof(fixedpoint_reltol_internal),typeof(resid),tTypeNoUnits,
                            typeof(integrator.tdir),typeof(integrator.k),typeof(sol),
                            typeof(integrator.rate_prototype),typeof(dde_f),
                            typeof(integrator.prog),typeof(integrator.cache),
                            typeof(integrator),typeof(prob),typeof(fixedpoint_norm),
                            typeof(integrator.opts),typeof(save_idxs)}(
                                sol, prob, integrator.u, integrator.k, integrator.t,
                                integrator.dt, dde_f, integrator.uprev, integrator.tprev,
                                u_cache, fixedpoint_abstol_internal,
                                fixedpoint_reltol_internal, resid, fixedpoint_norm,
                                alg.max_fixedpoint_iters, integrator.alg,
                                integrator.rate_prototype, integrator.notsaveat_idxs,
                                integrator.dtcache, integrator.dtchangeable,
                                integrator.dtpropose, integrator.tdir, integrator.EEst,
                                integrator.qold, integrator.q11, integrator.iter,
                                integrator.saveiter, integrator.saveiter_dense,
                                integrator.prog, integrator.cache, integrator.kshortsize,
                                integrator.just_hit_tstop, integrator.accept_step,
                                integrator.isout, integrator.reeval_fsal,
                                integrator.u_modified, integrator.opts, integrator,
                                save_idxs, saveat_internal, save_everystep, save_start)

    # set up additional initial values of newly created DDE integrator
    # (such as fsalfirst) and its callbacks
    initialize!(dde_int)
    initialize!(integrator.opts.callback, integrator.t, integrator.u, dde_int)

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

    # create vector of time points in the final solution
    n = integrator.save_everystep ? length(integrator.sol.t) : 2 # starting and end point
    n += length(integrator.saveat)
    integrator.save_start || (n -= 1)
    t_sol = Vector{typeof(integrator.t)}(n)

    # output initial time point if desired
    cur_ind = 1 # next index to write to
    if integrator.save_start
        t_sol[1] = integrator.sol.t[1]
        cur_ind = 2
    end

    # next indices to read from
    saveat_ind = 1
    sol_ind = 2 # already taken care of first element of solution of ODE integrator

    # merge additional time points and time points of solution of ODE integrator
    # assumed that both vectors are already sorted
    if integrator.save_everystep
        while saveat_ind <= length(integrator.saveat) && sol_ind <= length(integrator.sol.t)
            if integrator.tdir * integrator.saveat[saveat_ind] < integrator.tdir * integrator.sol.t[sol_ind]
                t_sol[cur_ind] = integrator.saveat[saveat_ind]
                saveat_ind += 1
            else
                t_sol[cur_ind] = integrator.sol.t[sol_ind]
                sol_ind += 1
            end

            cur_ind += 1
        end

        # copy remaining time points of solution of ODE integrator
        copy!(t_sol, cur_ind,
              integrator.sol.t, sol_ind, length(integrator.sol.t) + 1 - sol_ind)
    end
    # copy remaining additional time points
    copy!(t_sol, cur_ind,
          integrator.saveat, saveat_ind, length(integrator.saveat) + 1 - saveat_ind)

    # always output final time point
    t_sol[end] = integrator.sol.t[end]

    # evaluate solution at given time points
    u_sol = integrator.sol(t_sol; idxs=integrator.save_idxs).u

    # calculate analytical solutions to problem if existent
    if has_analytic(integrator.prob.f)
        if typeof(integrator.save_idxs) <: Void
            u_analytic = [integrator.prob.f(Val{:analytic}, t, integrator.sol[1])
                          for t in t_sol]
        else
            u_analytic = [integrator.prob.f(Val{:analytic}, t,
                                            integrator.sol[1])[integrator.save_idxs]
                          for t in t_sol]
        end
        errors = Dict{Symbol,recursive_eltype(u_sol)}()
    else
        u_analytic = nothing
        errors = nothing
    end

    # update interpolation data if only a subset of indices is returned
    if !(typeof(integrator.save_idxs) <: Void)
        if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
            interp = OrdinaryDiffEq.CompositeInterpolationData(
                integrator.sol.interp.f,
                [u[integrator.save_idxs] for u in integrator.sol.interp.timeseries],
                integrator.sol.interp.ts,
                [[k[integrator.save_idxs] for k in ks] for ks in integrator.sol.interp.ks],
                integrator.sol.interp.alg_choice, integrator.sol.interp.notsaveat_idxs,
                integrator.sol.interp.dense, integrator.sol.interp.cache)
        else
            interp = OrdinaryDiffEq.InterpolationData(
                integrator.sol.interp.f,
                [u[integrator.save_idxs] for u in integrator.sol.interp.timeseries],
                integrator.sol.interp.ts,
                [[k[integrator.save_idxs] for k in ks] for ks in integrator.sol.interp.ks],
                integrator.sol.interp.notsaveat_idxs, integrator.sol.interp.dense,
                integrator.sol.interp.cache)
        end
    else
        interp = integrator.sol.interp
    end

    # create updated solution which is evaluated at given time points and contains subset of
    # components and dense interpolation
    T = eltype(eltype(u_sol))
    N = length((size(u_sol[1])..., length(u_sol)))

    sol = ODESolution{T,N,typeof(u_sol),typeof(u_analytic),typeof(errors),typeof(t_sol),
                      typeof(interp.ks),typeof(integrator.sol.prob),
                      typeof(integrator.sol.alg),typeof(interp)}(
                          u_sol,u_analytic,errors,t_sol,interp.ks,integrator.sol.prob,
                          integrator.sol.alg,interp,integrator.sol.dense,
                          integrator.sol.tslocation,integrator.sol.retcode)

    # calculate errors of solution
    if sol.u_analytic != nothing
        calculate_solution_errors!(sol; fill_uanalytic=false,
                                   timeseries_errors=integrator.opts.timeseries_errors,
                                   dense_errors=integrator.opts.dense_errors)
    end

    sol.retcode = :Success
    return sol
end

function solve(prob::AbstractDDEProblem{uType,tType,lType,isinplace}, alg::algType,
               timeseries_init=uType[], ts_init=tType[], ks_init=[]; kwargs...) where
    {uType,tType,isinplace,algType<:AbstractMethodOfStepsAlgorithm,lType}

    integrator = init(prob, alg, timeseries_init, ts_init, ks_init; kwargs...)
    solve!(integrator)
end
