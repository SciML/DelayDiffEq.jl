using DelayDiffEq, OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = prob_dde_1delay_long

for constrained in (false, true)
    alg = MethodOfSteps(Tsit5(), constrained=constrained)
    sol = solve(prob, alg)

    @test allunique(sol.t)
end
