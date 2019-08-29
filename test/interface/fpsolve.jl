using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary, DiffEqDevTools
using LinearAlgebra
using Test

DDEProblemLibrary.importddeproblems()

const prob = DDEProblemLibrary.prob_dde_constant_2delays_long_ip
const prob_oop = DDEProblemLibrary.prob_dde_constant_2delays_long_oop
const prob_scalar = DDEProblemLibrary.prob_dde_constant_2delays_long_scalar

const testsol = TestSolution(solve(prob, MethodOfSteps(Vern9());
                             abstol = 1/10^14, reltol = 1/10^14))

@testset "NLFunctional" begin
  alg = MethodOfSteps(Tsit5(); fpsolve = NLFunctional(; max_iter = 10))

  ## in-place problem

  sol = solve(prob, alg)

  # check statistics
  @test sol.destats.nf > 3000
  @test sol.destats.nsolve == 0
  @test sol.destats.nfpiter > 300
  @test sol.destats.nfpconvfail > 50

  # compare it with the test solution
  sol2 = appxtrue(sol, testsol)
  @test sol2.errors[:L∞] < 2.7e-4

  ## out-of-place problem

  sol_oop = solve(prob_oop, alg)

  # compare it with the in-place solution
  @test sol_oop.destats.nf == sol.destats.nf
  @test sol_oop.destats.nsolve == sol.destats.nsolve
  @test sol_oop.destats.nfpiter == sol.destats.nfpiter
  @test sol_oop.destats.nfpconvfail == sol.destats.nfpconvfail
  @test sol_oop.t ≈ sol.t
  @test_broken sol_oop.u ≈ sol.u
  @test isapprox(sol.u, sol_oop.u; atol = 1e-7)

  ## scalar problem

  sol_scalar = solve(prob_scalar, alg)

  # compare it with the in-place solution
  @test sol_scalar.destats.nf == sol.destats.nf
  @test sol_scalar.destats.nsolve == sol.destats.nsolve
  @test sol_scalar.destats.nfpiter == sol.destats.nfpiter
  @test sol_scalar.destats.nfpconvfail == sol.destats.nfpconvfail
  @test sol_scalar.t ≈ sol.t
  @test sol_scalar.u ≈ sol[1, :]
end

@testset "NLAnderson" begin
  alg = MethodOfSteps(Tsit5(); fpsolve = NLAnderson(; max_iter = 10))

  ## in-place problem

  sol = solve(prob, alg)

  # check statistics
  @test sol.destats.nf < 2500
  @test sol.destats.nsolve > 0
  @test sol.destats.nfpiter < 250
  @test sol.destats.nfpconvfail < 50

  # compare it with the test solution
  sol2 = appxtrue(sol, testsol)
  @test sol2.errors[:L∞] < 2.7e-4

  ## out-of-place problem

  sol_oop = solve(prob_oop, alg)

  # compare it with the in-place solution
  @test_broken sol_oop.destats.nf == sol.destats.nf
  @test_broken sol_oop.destats.nsolve == sol.destats.nsolve
  @test_broken sol_oop.destats.nfpiter == sol.destats.nfpiter
  @test_broken sol_oop.destats.nfpconvfail == sol.destats.nfpconvfail
  @test_broken sol_oop.t ≈ sol.t
  @test_broken sol_oop.u ≈ sol.u
  @test appxtrue(sol, sol_oop).errors[:L∞] < 1.5e-7

  ## scalar problem
  
  sol_scalar = solve(prob_scalar, alg)

  # compare it with the in-place solution
  @test_broken sol_scalar.destats.nf == sol.destats.nf
  @test_broken sol_scalar.destats.nsolve == sol.destats.nsolve
  @test_broken sol_scalar.destats.nfpiter == sol.destats.nfpiter
  @test_broken sol_scalar.destats.nfpconvfail == sol.destats.nfpconvfail
  @test_broken sol_scalar.t ≈ sol.t
  @test_broken sol_scalar.u ≈ sol[1, :]
end