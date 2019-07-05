include("common.jl")

const prob_ip = prob_dde_constant_1delay_ip
const prob_scalar = prob_dde_constant_1delay_scalar
const ts = 0:0.1:10

# ODE algorithms
const algs = [Rosenbrock23(), Rosenbrock32(), ROS3P(), Rodas3(),
              RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
              Ros4LStab(), Rodas4(), Rodas42(), Rodas4P(), Rodas5()]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in algs
  println(nameof(typeof(alg)))

  stepsalg = MethodOfSteps(alg)
  sol_ip = solve(prob_ip, stepsalg)
  sol_scalar = solve(prob_scalar, stepsalg)

  @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)

  # Rodas4P tests fail on 32bit
  if Sys.WORD_SIZE == 32 && alg isa Rodas4P
    @test_broken sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
  else
    @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
  end
end
