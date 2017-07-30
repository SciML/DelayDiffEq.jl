function init(prob::AbstractDDEProblem{uType,tType,lType,isinplace}, alg::algType,
              timeseries_init=uType[], ts_init=tType[], ks_init=[];
              d_discontinuities=tType[], dtmax=tType(7*minimum(prob.lags)), dt=zero(tType),
              saveat=tType[], save_idxs=nothing, save_everystep=isempty(saveat),
              save_start=true, dense=save_everystep && !(typeof(alg) <: Discrete),
              minimal_solution=true, kwargs...) where
    {uType,tType,isinplace,algType<:AbstractMethodOfStepsAlgorithm,lType}

    # add lag locations to discontinuities vector
    d_discontinuities = [d_discontinuities; compute_discontinuity_tree(prob.lags, alg,
                                                                       prob.tspan[1])]

    # no fixed-point iterations for constrained algorithms,
    # and thus `dtmax` should match minimal lag
    if isconstrained(alg)
        dtmax = min(dtmax, prob.lags...)
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
    # assure that ODE integrator always satisfies tprev + dt == t
    integrator.dt = zero(integrator.dt)

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

    # create separate copies u and uprev, not pointing integrator.u or integrator.uprev,
    # containers for residuals and to cache uprev with correct dimensions and types
    # in particular for calculations with units residuals have to be unitless
    if typeof(integrator.u) <: AbstractArray
        u = recursivecopy(integrator.u)
        uprev = recursivecopy(integrator.uprev)
        uprev_cache = similar(integrator.u)
        resid = similar(integrator.u, uEltypeNoUnits)
    else
        u = deepcopy(integrator.u)
        uprev = deepcopy(integrator.uprev)
        resid = one(uEltypeNoUnits)
        uprev_cache = oneunit(eltype(uType))
    end

    # define caches for interpolation data
    k_cache = similar(integrator.k)
    k_integrator_cache = similar(integrator.k)

    # create heap of additional time points that will be contained in the solution
    # exclude the end point because of floating point issues and the starting point since it
    # is controlled by save_start
    if typeof(saveat) <: Number
        saveat_vec = collect(tType,
                             (prob.tspan[1] + integrator.tdir * abs(saveat)):
                             integrator.tdir * abs(saveat):
                             (prob.tspan[end] - integrator.tdir * abs(saveat)))
    else
        saveat_vec = collect(tType, Iterators.filter(
            x -> integrator.tdir * prob.tspan[1] < integrator.tdir * x < integrator.tdir *
            prob.tspan[end],
            saveat))
    end

    if integrator.tdir > 0
        saveat_internal = binary_minheap(saveat_vec)
    else
        saveat_internal = binary_maxheap(saveat_vec)
    end

    # check if all indices should be returned
    if !(typeof(save_idxs) <: Void) && collect(save_idxs) == collect(1:length(integrator.u))
        save_idxs = nothing # prevents indexing of ODE solution and saves memory
    end

    # separate options of integrator and ODE integrator since ODE integrator always saves
    # every step and every index (necessary for history function)
    opts = OrdinaryDiffEq.DEOptions(integrator.opts.maxiters,
                                    integrator.opts.timeseries_steps, save_everystep,
                                    integrator.opts.adaptive, integrator.opts.abstol,
                                    integrator.opts.reltol,
                                    tTypeNoUnits(integrator.opts.gamma),
                                    tTypeNoUnits(integrator.opts.qmax),
                                    tTypeNoUnits(integrator.opts.qmin),
                                    tType(integrator.opts.dtmax),
                                    tType(integrator.opts.dtmin),
                                    integrator.opts.internalnorm, save_idxs,
                                    integrator.opts.tstops, saveat_internal,
                                    integrator.opts.d_discontinuities,
                                    integrator.opts.userdata, integrator.opts.progress,
                                    integrator.opts.progress_steps,
                                    integrator.opts.progress_name,
                                    integrator.opts.progress_message,
                                    integrator.opts.timeseries_errors,
                                    integrator.opts.dense_errors,
                                    tTypeNoUnits(integrator.opts.beta1),
                                    tTypeNoUnits(integrator.opts.beta2),
                                    tTypeNoUnits(integrator.opts.qoldinit),
                                    dense && integrator.opts.dense, save_start,
                                    integrator.opts.callback, integrator.opts.isoutofdomain,
                                    integrator.opts.unstable_check,
                                    integrator.opts.verbose, integrator.opts.calck,
                                    integrator.opts.force_dtmin,
                                    integrator.opts.advance_to_tstop,
                                    integrator.opts.stop_at_next_tstop)

    # need copy of heap of additional time points (nodes will be deleted!) in order to
    # remove unneeded time points of ODE solution as soon as possible and keep track
    # of passed time points
    saveat_copy = deepcopy(opts.saveat)

    # reduction of solution only possible if no dense interpolation required and only
    # selected time points saved
    minimal_solution = minimal_solution && !opts.dense && !opts.save_everystep

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
                            typeof(opts),typeof(saveat_copy)}(
                                sol, prob, u, integrator.k, integrator.t, dt, dde_f, uprev,
                                integrator.tprev, uprev_cache, k_cache, k_integrator_cache,
                                fixedpoint_abstol_internal, fixedpoint_reltol_internal,
                                resid, fixedpoint_norm, alg.max_fixedpoint_iters,
                                minimal_solution, integrator.alg, integrator.rate_prototype,
                                integrator.notsaveat_idxs, integrator.dtcache,
                                integrator.dtchangeable, integrator.dtpropose,
                                integrator.tdir, integrator.EEst, integrator.qold,
                                integrator.q11, integrator.iter, integrator.saveiter,
                                integrator.saveiter_dense, integrator.prog,
                                integrator.cache, integrator.kshortsize,
                                integrator.just_hit_tstop, integrator.accept_step,
                                integrator.isout, integrator.reeval_fsal,
                                integrator.u_modified, opts, integrator, saveat_copy)

    # set up additional initial values of newly created DDE integrator
    # (such as fsalfirst) and its callbacks
    initialize!(dde_int)
    initialize!(integrator.opts.callback, integrator.t, u, dde_int)

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

    # create array of time points and values that form solution
    sol_array = build_solution_array(integrator)

    # create interpolation data of solution
    interp = build_solution_interpolation(integrator, sol_array)

    # calculate analytical solutions to problem if existent
    if has_analytic(integrator.prob.f)
        if typeof(integrator.opts.save_idxs) <: Void
            u_analytic = [integrator.prob.f(Val{:analytic}, t, integrator.sol[1])
                          for t in sol_array.t]
        else
            u_analytic = [@view(integrator.prob.f(
                Val{:analytic}, t, integrator.sol[1])[integrator.opts.save_idxs])
                          for t in sol_array.t]
        end
        errors = Dict{Symbol,eltype(integrator.u)}()
    else
        u_analytic = nothing
        errors = nothing
    end

    # combine arrays of time points and values, interpolation data, and analytical solution
    # to solution
    sol = ODESolution{eltype(sol_array),size(sol_array),typeof(sol_array.u),
                      typeof(u_analytic),typeof(errors),typeof(sol_array.t),
                      typeof(interp.ks),typeof(integrator.sol.prob),
                      typeof(integrator.sol.alg),typeof(interp)}(
                          sol_array.u, u_analytic, errors, sol_array.t, interp.ks,
                          integrator.sol.prob, integrator.sol.alg, interp, interp.dense,
                          integrator.sol.tslocation, integrator.sol.retcode)

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
