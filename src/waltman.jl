using CSV
using DataFrames
using DelayDiffEq
using DDEProblemLibrary

"""
    waltman()

Solve the stiff delay-differential equation problem `DDEProblemLibrary.prob_dde_RADAR5_waltman_5` for time points `0:0.1:300`.

The simulation is computed with `MethodOfSteps(KenCarp5())` and the tolerances used by Guglielmi and Hairer (Figure 3): the relative tolerance is set to `1e-9`, and absolute tolerances are `1e-21` for the first four components and `1e-9` for the fifth and sixth component.

## References

Guglielmi, N. and E. Hairer. (2001). Implementing Radau IIA methods for stiff delay differential equations, Computing 67(1), pp. 1-12.

Waltman, P. (1978). A threshold model of antigen-stimulated antibody production, Theoretical Immunology (8), pp. 437-453.
"""
function waltman()
    # DDE problem and tolerances used by Guglielmi and Hairer (Figure 3)
    prob = DDEProblemLibrary.prob_dde_RADAR5_waltman_5
    reltol = 1e-9
    abstol = [
        1e-21,
        1e-21,
        1e-21,
        1e-21,
        1e-9,
        1e-9,
    ]

    # Compute solution
    alg = MethodOfSteps(KenCarp5())
    sol = solve(prob, alg; reltol = reltol,
                abstol = abstol,
                saveat = 0.1)

    return sol
end

"""
    waltman(io)

Write the solution computed by [`waltman()`](@ref) to `io` in the CSV file format, and return the solution object.

!!! note
    If `io` is a file name, its directory is created if it does not exist yet.
"""
function waltman(io)
    # Solve the DDE problem
    sol = waltman()

    # Create a dataframe of the solution
    df = DataFrame(sol)

    # Output the dataframe in the CSV file format
    io isa AbstractString && mkpath(dirname(io))
    CSV.write(io, df)

    return sol
end
