include("common.jl")

@testset for constrained in (false, true)
  sol = solve(prob_dde_constant_1delay_long_ip,
              MethodOfSteps(Tsit5(), constrained=constrained))

  @test allunique(sol.t)
end
