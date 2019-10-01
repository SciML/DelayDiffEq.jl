using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using LinearAlgebra, Test

DDEProblemLibrary.importddeproblems()
const PROB_WALTMAN = DDEProblemLibrary.prob_dde_RADAR5_waltman_5
const PROB_KWARGS = (reltol = 1e-6, abstol = [1e-21, 1e-21, 1e-21, 1e-21, 1e-9, 1e-9])

# solution at final time point T = 300 obtained from RADAR5
# with relative tolerance 1e-6 and componentwise absolute tolerances
# 1e-21, 1e-21, 1e-21, 1e-21, 1e-9, and 1e-9
const RADAR5_SOL = [6.154488183e-16, 3.377120916e-7, 4.22140331e-7,
                    2.142554562e-6, 299.9999999, 299.6430338]

function test_waltman_sol(sol)
    @test sol.retcode === :Success
    @test sol.t[end] == 300

    # compare solution at the final time point with RADAR5
    for i in 1:6
        @test sol.u[end][i] ≈ RADAR5_SOL[i] rtol = 1e-3 atol = 1e-17
    end
end

# standard factorization
sol1 = solve(PROB_WALTMAN, MethodOfSteps(Rosenbrock23()); PROB_KWARGS...)
test_waltman_sol(sol1)

# in-place LU factorization
sol2 = solve(PROB_WALTMAN, MethodOfSteps(Rosenbrock23(linsolve=LinSolveFactorize(lu!)));
             PROB_KWARGS...)
test_waltman_sol(sol2)

# out-of-place LU factorization
sol3 = solve(PROB_WALTMAN, MethodOfSteps(Rosenbrock23(linsolve=LinSolveFactorize(lu)));
             PROB_KWARGS...)
test_waltman_sol(sol3)

# compare in-place and out-of-place LU factorization
@test sol2.t == sol3.t
@test sol2.u == sol3.u
