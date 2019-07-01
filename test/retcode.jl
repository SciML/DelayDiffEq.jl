include("common.jl")

@testset for composite in (true, false)
  alg = MethodOfSteps(composite ? AutoTsit5(Rosenbrock23()) : Tsit5();
                      constrained=false)

  sol1 = solve(prob_dde_constant_1delay_ip, alg)
  @test sol1.retcode == :Success

  sol2 = solve(prob_dde_constant_1delay_ip, alg; maxiters = 1, verbose = false)
  @test sol2.retcode == :MaxIters

  sol3 = solve(prob_dde_constant_1delay_ip, alg; dtmin = 5, verbose = false)
  @test sol3.retcode == :DtLessThanMin
end
