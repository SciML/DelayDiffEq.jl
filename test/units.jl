include("common.jl")
using DiffEqProblemLibrary.DDEProblemLibrary: remake_dde_constant_u0_tType
using Unitful

const probs =
  Dict(true => remake_dde_constant_u0_tType(prob_dde_constant_1delay_long_ip, [1.0u"N"],
                                            typeof(1.0u"s")),
       false => remake_dde_constant_u0_tType(prob_dde_constant_1delay_long_scalar, 1.0u"N",
                                             typeof(1.0u"s")))

# we test the current handling of units for regressions
# however, it is broken upstream: https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/issues/828

@testset for inplace in (true, false)
  prob = probs[inplace]

  # default tolerances
  alg1 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100)
  sol1 = solve(prob, alg1)

  # without units
  alg2 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                       fixedpoint_abstol=1e-6, fixedpoint_reltol=1e-3)

  if inplace
    @test_throws Unitful.DimensionError solve(prob, alg2)
  else
    sol2 = solve(prob, alg2)

    @test sol1.t == sol2.t
    @test sol1.u == sol2.u
  end

  # with correct units
  alg3 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                       fixedpoint_abstol=1e-6u"N", fixedpoint_reltol=1e-3u"N")
  sol3 = solve(prob, alg3)

  @test sol1.t == sol3.t
  @test sol1.u == sol3.u

  # with correct units as vectors
  if inplace
    alg4 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                         fixedpoint_abstol=[1e-6u"N"], fixedpoint_reltol=[1e-3u"N"])
    sol4 = solve(prob, alg4)

    @test sol1.t == sol4.t
    @test sol1.u == sol4.u
  end

  # with incorrect units for absolute tolerance
  alg5 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                       fixedpoint_abstol=1e-6u"s", fixedpoint_reltol=1e-3u"N")
  @test_throws Unitful.DimensionError solve(prob, alg5)

  # with incorrect units for relative tolerance
  alg6 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                       fixedpoint_abstol=1e-6u"N", fixedpoint_reltol=1e-3u"s")
  @test_throws Unitful.DimensionError solve(prob, alg6)
end
