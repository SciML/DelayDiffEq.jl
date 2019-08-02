using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using Test

DDEProblemLibrary.importddeproblems()

const prob = DDEProblemLibrary.prob_dde_constant_1delay_ip

@testset for composite in (true, false)
  alg = MethodOfSteps(composite ? AutoTsit5(Rosenbrock23()) : Tsit5();
                      constrained=false)

  sol1 = solve(prob, alg)
  @test sol1.retcode == :Success

  sol2 = solve(prob, alg; maxiters = 1, verbose = false)
  @test sol2.retcode == :MaxIters

  sol3 = solve(prob, alg; dtmin = 5, verbose = false)
  @test sol3.retcode == :DtLessThanMin
end
