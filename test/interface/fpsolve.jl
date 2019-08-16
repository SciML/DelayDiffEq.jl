using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using Test

DDEProblemLibrary.importddeproblems()

const prob = DDEProblemLibrary.prob_dde_constant_2delays_ip
const prob_oop = DDEProblemLibrary.prob_dde_constant_2delays_oop
const prob_scalar = DDEProblemLibrary.prob_dde_constant_2delays_scalar

@testset "NLFunctional" begin
  alg = MethodOfSteps(Tsit5(); fpsolve = NLFunctional(; max_iter = 10))

  sol = solve(prob, alg)
  @test sol.errors[:l∞] < 4.5e-3

  sol_oop = solve(prob_oop, alg)
  @test sol.t ≈ sol_oop.t
  @test sol.u ≈ sol_oop.u

  sol_scalar = solve(prob_scalar, alg)
  @test sol.t ≈ sol_scalar.t
  @test sol[1, :] ≈ sol_scalar.u
end

@testset "NLAnderson" begin
  alg = MethodOfSteps(Tsit5(); fpsolve = NLAnderson(; max_iter = 10))

  sol = solve(prob, alg)
  @test sol.errors[:l∞] < 4.5e-3

  sol_oop = solve(prob_oop, alg)
  @test sol.t ≈ sol_oop.t
  @test sol.u ≈ sol_oop.u

  sol_scalar = solve(prob_scalar, alg)
  @test sol.t ≈ sol_scalar.t
  @test sol[1, :] ≈ sol_scalar.u
end