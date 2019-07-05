include("common.jl")

const prob_ip = prob_dde_constant_1delay_ip
const prob_scalar = prob_dde_constant_1delay_scalar
const ts = 0:0.1:10

# ODE algorithms
const algs = [GenericImplicitEuler(), GenericTrapezoid(),
              ImplicitEuler(), ImplicitMidpoint(), Trapezoid(),
              TRBDF2(), SDIRK2(), SSPSDIRK2(),
              Kvaerno3(), KenCarp3(),
              Cash4(), Hairer4(), Hairer42(), Kvaerno4(), KenCarp4(),
              Kvaerno5(), KenCarp5()]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in algs
  println(nameof(typeof(alg)))

  stepsalg = MethodOfSteps(alg)
  sol_ip = solve(prob_ip, stepsalg)
  sol_scalar = solve(prob_scalar, stepsalg)

  @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
  @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
end
