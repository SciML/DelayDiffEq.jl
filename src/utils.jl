"""
    compute_discontinuity_tree(lags, alg, start_val)

Compute time points of discontinuities for delays `lags` up to the order of algorithm `alg`,
starting at initial time point `start_val`.

# Examples

```jldoctest
julia> compute_discontinuity_tree([1//2], BS3(), 1)
3-element Array{Rational{Int64},1}:
 3//2
 2//1
 5//2

julia> compute_discontinuity_tree([1//2, 1//3], BS3(), 1)
8-element Array{Rational{Int64},1}:
  3//2
  4//3
  2//1
 11//6
  5//3
  5//2
  7//3
 13//6
```
"""
function compute_discontinuity_tree(lags, alg, start_val)
    start_val + unique(vcat((sum.(collect(with_replacement_combinations(lags, i))) for i in
                             1:alg_order(alg))...))
end

"""
    reduce_solution!(integrator::DDEIntegrator, tmax)

Remove time points of ODE solution of `integrator` up to time `tmax` that are not required
for calculation of DDE solution.
"""
@inline function reduce_solution!(integrator::DDEIntegrator, tmax)
    if !integrator.opts.dense && !integrator.opts.save_everystep &&
        integrator.saveiter < integrator.integrator.saveiter

        # last time point of ODE solution that is known to be required for DDE solution
        @inbounds t_sol_prev = integrator.integrator.sol.t[integrator.saveiter]
        # next time point of ODE solution that might be required for DDE solution
        @inbounds t_sol = integrator.integrator.sol.t[integrator.saveiter + 1]

        @inbounds while t_sol < tmax && integrator.saveiter < integrator.integrator.saveiter
            needed = false

            while !(isempty(integrator.saveat)) && top(integrator.saveat) <= t_sol
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
            needed = needed || (!(isempty(integrator.saveat)) &&
                                integrator.saveiter + 1 < integrator.integrator.saveiter &&
                                t_sol < top(integrator.saveat) <
                                integrator.integrator.sol.t[integrator.saveiter + 2])

            if !needed
                # delete not required time point, function value, and interpolation data
                deleteat!(integrator.integrator.sol.t, integrator.saveiter + 1)
                deleteat!(integrator.integrator.sol.u, integrator.saveiter + 1)
                deleteat!(integrator.integrator.sol.k, integrator.saveiter + 1)

                # correct other arrays
                pop!(integrator.integrator.notsaveat_idxs)
                if typeof(integrator.integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
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
end
