using DelayDiffEq, Test

@testset "in-place" begin
  # define functions (Hutchinson's equation)
  function f(du, u, h, p, t)
    du[1] = u[1] * (1  - h(p, t - 1)[1])
    nothing
  end

  function g(J, u, h, p, t)
    J[1, 1] = 1 - h(p, t - 1)[1]
    nothing
  end

  h(p, t) = [0.0]

  # define problems
  prob_wo_jac = DDEProblem(DDEFunction{true}(f), [1.0], h, (0.0, 40.0);
                           constant_lags = [1])
  prob_w_jac = DDEProblem(DDEFunction{true}(f; jac = g), [1.0], h, (0.0, 40.0);
                          constant_lags = [1])

  # compute solutions
  for alg in (Rosenbrock23(), TRBDF2())
    sol_wo_jac = solve(prob_wo_jac, MethodOfSteps(alg))
    sol_w_jac = solve(prob_w_jac, MethodOfSteps(alg))

    @test sol_wo_jac.t ≈ sol_w_jac.t
    @test sol_wo_jac.u ≈ sol_w_jac.u
  end
end

@testset "out-of-place" begin
  # define functions (Hutchinson's equation)
  f(u, h, p, t) = [u[1] * (1  - h(p, t - 1)[1])]

  g(u, h, p, t) = fill(1 - h(p, t - 1)[1], 1, 1)

  h(p, t) = [0.0]

  # define problems
  prob_wo_jac = DDEProblem(DDEFunction{false}(f), [1.0], h, (0.0, 40.0);
                           constant_lags = [1])
  prob_w_jac = DDEProblem(DDEFunction{false}(f; jac = g), [1.0], h, (0.0, 40.0);
                          constant_lags = [1])

  # compute solutions
  for alg in (Rosenbrock23(), TRBDF2())
    sol_wo_jac = solve(prob_wo_jac, MethodOfSteps(alg))
    sol_w_jac = solve(prob_w_jac, MethodOfSteps(alg))

    @test sol_wo_jac.t ≈ sol_w_jac.t
    @test sol_wo_jac.u ≈ sol_w_jac.u
  end
end
