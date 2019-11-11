using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using Test

DDEProblemLibrary.importddeproblems()

const prob_ip = DDEProblemLibrary.prob_dde_constant_1delay_ip
const prob_scalar = DDEProblemLibrary.prob_dde_constant_1delay_scalar
const ts = 0:0.1:10

# ODE algorithms
const algs = Dict(ImplicitEuler() => false,
                  ImplicitMidpoint() => true,
                  Trapezoid() => false,
                  TRBDF2() => false,
                  SDIRK2() => false,
                  SSPSDIRK2() => true,
                  Kvaerno3() => false,
                  KenCarp3() => false,
                  Cash4() => false,
                  Hairer4() => false,
                  Hairer42() => false,
                  Kvaerno4() => false,
                  KenCarp4() => false,
                  Kvaerno5() => false,
                  KenCarp5() => true)

@testset "Algorithm $(nameof(typeof(alg)))" for (alg, pass) in algs
  println(nameof(typeof(alg)))

  stepsalg = MethodOfSteps(alg)
  sol_ip = solve(prob_ip, stepsalg; dt = 0.1)
  sol_scalar = solve(prob_scalar, stepsalg; dt = 0.1)

  if pass
    @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
    @test sol_ip.t ≈ sol_scalar.t
    @test sol_ip[1, :] ≈ sol_scalar.u
  else
    @test_broken sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
    @test_broken sol_ip.t ≈ sol_scalar.t
    @test_broken sol_ip[1, :] ≈ sol_scalar.u
  end
end
