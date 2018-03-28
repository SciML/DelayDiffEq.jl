@testset "Return code" begin
    alg = MethodOfSteps(BS3(); constrained=false)

    sol1 = solve(prob_dde_1delay, alg)
    @test sol1.retcode == :Success

    sol2 = solve(prob_dde_1delay, alg; maxiters = 1, verbose = false)
    @test sol2.retcode == :MaxIters

    sol3 = solve(prob_dde_1delay, alg; dtmin = 5, verbose = false)
    @test sol3.retcode == :DtLessThanMin
end
