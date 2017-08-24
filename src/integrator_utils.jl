"""
    reduce_solution!(integrator::DDEIntegrator, tmax)

Remove time points of ODE solution of `integrator` up to time `tmax` that are not required
for calculation of DDE solution.
"""
@inline function reduce_solution!(integrator::DDEIntegrator, tmax)
    if integrator.minimal_solution
        if integrator.saveiter < integrator.integrator.saveiter
            # last time point of ODE solution that is known to be required for DDE solution
            @inbounds t_sol_prev = integrator.integrator.sol.t[integrator.saveiter]
            # next time point of ODE solution that might be required for DDE solution
            @inbounds t_sol = integrator.integrator.sol.t[integrator.saveiter + 1]

            @inbounds while t_sol < tmax && integrator.saveiter <
                integrator.integrator.saveiter

                needed = false

                while !isempty(integrator.saveat) && top(integrator.saveat) <= t_sol
                    # do not remove time point if it is upper bound of a required
                    # interpolation interval
                    if t_sol_prev < top(integrator.saveat)
                        needed = true
                    end

                    # remove passed time points
                    pop!(integrator.saveat)
                end

                # do not remove time point if it is lower bound of a required interpolation
                # interval
                needed = needed || (!isempty(integrator.saveat) &&
                                    integrator.saveiter + 1 <
                                    integrator.integrator.saveiter &&
                                    t_sol < top(integrator.saveat) <
                                    integrator.integrator.sol.t[integrator.saveiter + 2])

                if !needed
                    # delete not required time point, function value, and interpolation data
                    deleteat!(integrator.integrator.sol.t, integrator.saveiter + 1)
                    deleteat!(integrator.integrator.sol.u, integrator.saveiter + 1)
                    deleteat!(integrator.integrator.sol.k, integrator.saveiter + 1)

                    # correct other arrays
                    pop!(integrator.integrator.notsaveat_idxs)
                    if typeof(integrator.integrator.alg) <: OrdinaryDiffEq.OrdinaryDiffEqCompositeAlgorithm
                        deleteat!(integrator.integrator.alg_choice, integrator.saveiter + 1)
                    end

                    # update counter of saved time points
                    integrator.integrator.saveiter -= 1
                    integrator.integrator.saveiter_dense -= 1

                    # update time point of loop
                    t_sol = integrator.integrator.sol.t[integrator.saveiter + 1]
                else
                    # increase counter of time points that are required for DDE solution
                    integrator.saveiter += 1

                    # update time points of loop
                    t_sol_prev = t_sol
                    t_sol = integrator.integrator.sol.t[integrator.saveiter + 1]
                end
            end
        end
    else
        # do not reduce ODE solution and only remove passed time points up to time tmax
        # this ensures that no time points after tmax are added to the DDE solution
        while !isempty(integrator.saveat) && top(integrator.saveat) <= tmax
            pop!(integrator.saveat)
        end
    end
end

"""
    build_solution_array(integrator::DDEIntegrator)

Create a `DiffEqArray` of the time points and values that form the solution of `integrator`.
"""
function build_solution_array(integrator::DDEIntegrator)
    if integrator.opts.save_everystep && isempty(integrator.opts.saveat)
        # use solution of ODE integrator if no additional time points provided
        if integrator.opts.save_start
            t = integrator.sol.t
            if typeof(integrator.opts.save_idxs) <: Void
                u = integrator.sol.u
            else
                u = [@view(u[integrator.opts.save_idxs]) for u in integrator.sol.u]
            end
        else # remove initial time point
            t = @view integrator.sol.t[2:end]
            if typeof(integrator.opts.save_idxs) <: Void
                u = @view integrator.sol.u[2:end]
            else
                u = [@view(u[integrator.opts.save_idxs]) for u in
                     Iterators.drop(integrator.sol.u, 1)]
            end
        end
    else
        # calculate number of additional time points of solution
        saveat_length = length(integrator.opts.saveat) - length(integrator.saveat)

        # create vectors of time points and corresponding values which form final solution
        n = integrator.opts.save_everystep ? length(integrator.sol.t) : 2
        n += saveat_length
        integrator.opts.save_start || (n -= 1)
        t = Vector{typeof(integrator.t)}(n)
        u = Vector{typeof(integrator.u)}(n)

        # output initial time point if desired
        write_idx = 1 # next index of solution to write to
        if integrator.opts.save_start
            t[1] = integrator.sol.t[1]
            if typeof(integrator.opts.save_idxs) <: Void
                u[1] = integrator.sol.u[1]
            else
                u[1] = @view integrator.sol.u[1][integrator.opts.save_idxs]
            end
            write_idx = 2
        end

        # merge additional time points and time points of solution of ODE integrator
        # both data structures are already sorted
        if integrator.opts.save_everystep
            sol_idx = 2 # next index of solution of ODE integrator to read from
            sol_length = length(integrator.sol.t)

            @inbounds while saveat_length > 0 && sol_idx < sol_length
                if integrator.tdir * top(integrator.opts.saveat) <
                    integrator.tdir * integrator.sol.t[sol_idx]

                    # copy additional time points and calculate corresponding values
                    t[write_idx] = pop!(integrator.opts.saveat)
                    u[write_idx] = integrator.sol(t[write_idx];
                                                  idxs=integrator.opts.save_idxs)

                    saveat_length -= 1
                else
                    # copy time points and values of solution of ODE integrator
                    t[write_idx] = integrator.sol.t[sol_idx]
                    if typeof(integrator.opts.save_idxs) <: Void
                        u[write_idx] = integrator.sol.u[sol_idx]
                    else
                        u[write_idx] =
                            @view integrator.sol.u[sol_idx][integrator.opts.save_idxs]
                    end

                    sol_idx += 1
                end

                write_idx += 1
            end

            # copy remaining time points of solution of ODE integrator
            # except of final time point
            copy!(t, write_idx, integrator.sol.t, sol_idx, sol_length - sol_idx)
            if typeof(integrator.opts.save_idxs) <: Void
                copy!(u, write_idx, integrator.sol.u, sol_idx, sol_length - sol_idx)
            else
                copy!(u, write_idx,
                      (@view(u[integrator.opts.save_indxs]) for u in integrator.sol.u),
                      sol_idx, sol_length - sol_idx)
            end

            write_idx += sol_length - sol_idx
        end

        # copy remaining additional time points and calculate corresponding values
        @inbounds while saveat_length > 0
            t[write_idx] = pop!(integrator.opts.saveat)
            u[write_idx] = integrator.sol(t[write_idx]; idxs=integrator.opts.save_idxs)

            saveat_length -= 1
            write_idx += 1
        end

        # always output final time point
        t[end] = integrator.sol.t[end]
        if typeof(integrator.opts.save_idxs) <: Void
            u[end] = integrator.sol.u[end]
        else
            u[end] = @view integrator.sol.u[end][integrator.opts.save_idxs]
        end
    end

    DiffEqArray(u, t)
end

"""
    build_solution_interpolation(integrator::DDEIntegrator, sol::DiffEqArray)

Create interpolation data to solution of `integrator`, which is formed by time points and
values in `sol`.
"""
function build_solution_interpolation(integrator::DDEIntegrator, sol::DiffEqArray)
    if integrator.opts.dense
        if typeof(integrator.opts.save_idxs) <: Void
            integrator.sol.interp
        else # update interpolation data if only a subset of indices is returned
            if typeof(integrator.alg) <: OrdinaryDiffEq.OrdinaryDiffEqCompositeAlgorithm
                OrdinaryDiffEq.CompositeInterpolationData(
                    integrator.sol.interp.f, [@view(u[integrator.opts.save_idxs]) for u in
                                              integrator.sol.interp.timeseries],
                    integrator.sol.interp.ts, [[@view(k[integrator.opts.save_idxs]) for k in
                                                ks] for ks in integrator.sol.interp.ks],
                    integrator.sol.interp.alg_choice, integrator.sol.interp.notsaveat_idxs,
                    true, integrator.sol.interp.cache)
            else
                OrdinaryDiffEq.InterpolationData(
                    integrator.sol.interp.f, [@view(u[integrator.opts.save_idxs]) for u in
                                              integrator.sol.interp.timeseries],
                    integrator.sol.interp.ts, [[@view(k[integrator.opts.save_idxs]) for k in
                                                ks] for ks in integrator.sol.interp.ks],
                    integrator.sol.interp.notsaveat_idxs, true,
                    integrator.sol.interp.cache)
            end
        end
    else # create not dense interpolation data if desired
        if typeof(integrator.alg) <: OrdinaryDiffEq.OrdinaryDiffEqCompositeAlgorithm
            OrdinaryDiffEq.CompositeInterpolationData(
                integrator.sol.interp.f, sol.u, sol.t, typeof(integrator.sol.k)(0), Int[],
                Int[], false, integrator.sol.interp.cache)
        else
            OrdinaryDiffEq.InterpolationData(
                integrator.sol.interp.f, sol.u, sol.t, typeof(integrator.sol.k)(0), Int[],
                false, integrator.sol.interp.cache)
        end
    end
end

"""
    update_ode_integrator!(integrator::DDEIntegrator)

Update ODE integrator of `integrator` to current time interval, values and interpolation
data of `integrator`.
"""
function update_ode_integrator!(integrator::DDEIntegrator)
    # update time interval of ODE integrator
    integrator.integrator.t = integrator.t
    integrator.integrator.tprev = integrator.tprev
    integrator.integrator.dt = integrator.dt

    # copy u(tprev) since it is overwritten by integrator at the end of apply_step!
    if typeof(integrator.u) <: AbstractArray
        recursivecopy!(integrator.integrator.u, integrator.u)
        recursivecopy!(integrator.integrator.uprev, integrator.uprev)
    else
        integrator.integrator.u = integrator.u
        integrator.integrator.uprev = integrator.uprev
    end

    # add additional interpolation steps
    OrdinaryDiffEq.ode_addsteps!(integrator)

    # copy interpolation data (fsalfirst overwritten at the end of apply_step!, which also
    # updates k[1] when using chaches for which k[1] points to fsalfirst)
    recursivecopy!(integrator.integrator.k, view(integrator.k, 1:integrator.kshortsize))
    @tight_loop_macros for i in integrator.kshortsize+1:length(integrator.k)
        @inbounds copyat_or_push!(integrator.integrator.k, i, integrator.k[i], Val{false})
    end

    # remove additional interpolation steps
    resize!(integrator.k, integrator.kshortsize)
end
