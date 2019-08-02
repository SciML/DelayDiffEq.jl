using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using Test

DDEProblemLibrary.importddeproblems()

const prob = DDEProblemLibrary.prob_dde_constant_1delay_long_ip

@testset for constrained in (false, true)
  sol = solve(prob, MethodOfSteps(Tsit5(); constrained = constrained))

  @test allunique(sol.t)
end
