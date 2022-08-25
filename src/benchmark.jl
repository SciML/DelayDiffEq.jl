using CSV
using DataFrames
using DelayDiffEq
using DiffEqDevTools
using DDEProblemLibrary

function benchmark()
    # Compute reference solution of third state
    # with low tolerances
    fpsolve = NLFunctional(; max_iter = 1000)
    alg = MethodOfSteps(Vern9(); fpsolve)
    sol = TestSolution(solve(prob_dde_qs, alg;
                             reltol = 1e-14,
                             abstol = 1e-14,
                             save_idxs = 3))

    # Benchmark different algorithms
    setups = [
        Dict(:alg => MethodOfSteps(RK4())),
        Dict(:alg => MethodOfSteps(DP5())),
        Dict(:alg => MethodOfSteps(OwrenZen3())),
        Dict(:alg => MethodOfSteps(OwrenZen5())),
        Dict(:alg => MethodOfSteps(Vern6())),
        Dict(:alg => MethodOfSteps(Vern7())),
        Dict(:alg => MethodOfSteps(Vern8())),
        Dict(:alg => MethodOfSteps(Vern9())),
        Dict(:alg => MethodOfSteps(Rosenbrock23())),
        Dict(:alg => MethodOfSteps(Rodas4())),
        Dict(:alg => MethodOfSteps(Rodas5P())),
    ]
    abstols = exp10.(-(4:11))
    reltols = exp10.(-(1:8))
    return WorkPrecisionSet(prob_dde_qs,
                            abstols,
                            reltols,
                            setups;
                            save_idxs = 3,
                            appxsol = sol,
                            maxiters = Int(1e5),
                            error_estimate = :l2)
end

"""
    benchmark(io)

Write the results of [`benchmark()`](@ref) to `io` in the CSV file format, and return the results.

!!! note
    If `io` is a file name, its directory is created if it does not exist yet.
"""
function benchmark(io)
    # Run benchmark
    wp = benchmark()

    # Create long data frame 
    df = DataFrame(; alg = String[],
                   error = Float64[],
                   time = Float64[])
    for (i, alg) in enumerate(wp.names)
        wp_i = wp[i]
        append!(df,
                DataFrame(; alg = alg,
                          error = wp_i.errors,
                          time = wp_i.times))
    end

    # Write data frame in CSV file format
    io isa AbstractString && mkpath(dirname(io))
    CSV.write(io, df)

    return wp
end

## Separate documentation
## (fixes issues with listing package in LaTeX)

@doc """
    benchmark()

Benchmark the delay-differential equation problem
`DDEProblemLibrary.prob_dde_qs` of quorum sensing in
/Pseudomonas putida/.

The benchmarked algorithms are:
- `MethodOfSteps(RK4())`
- `MethodOfSteps(DP5())`
- `MethodOfSteps(OwrenZen3())`
- `MethodOfSteps(OwrenZen5())`
- `MethodOfSteps(Vern6())`
- `MethodOfSteps(Vern7())`
- `MethodOfSteps(Vern8())`
- `MethodOfSteps(Vern9())`
- `MethodOfSteps(Rosenbrock23())`
- `MethodOfSteps(Rodas4())`
- `MethodOfSteps(Rodas5P())`

## Reference

Buddrus-Schiemann, K. et al. (2014) “Analysis of N-acylhomoserine lactone dynamics in continuous cultures of Pseudomonas putida IsoF by use of ELISA and UHPLC/qTOF-MS-derived measurements and mathematical models,” Analytical and Bioanalytical Chemistry. Springer Science and Business Media LLC. doi: [10.1007/s00216-014-8063-6](https://doi.org/10.1007/s00216-014-8063-6). 
""" benchmark_qs
