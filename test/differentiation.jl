using DelayDiffEq, Calculus, ForwardDiff
using LinearAlgebra, Test

# Hutchinson's equation
function f(du, u, h, p, t)
  du[1] = p[1] * u[1] * (1 - h(p, t - p[2])[1] / p[3])
  nothing
end

h(p, t) = [p[4]]

@testset "Gradient" begin
  function test(p)
    prob = DDEProblem(f, [p[5]], h, eltype(p).((0.0, 20.0)), copy(p);
                      constant_lags = [p[2]])
    sol = solve(prob, MethodOfSteps(Tsit5()); abstol = 1e-14, reltol = 1e-14)
    diff = @. sol[1, :] - 10 * exp(-sol.t)
    dot(diff, diff)
  end

  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = Calculus.finite_difference(test, p)
  fordiff = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.gradient(test, p)

  @test_broken findiff ≈ fordiff
end

@testset "Jacobian" begin
  function test(p)
    prob = DDEProblem(f, [p[5]], h, eltype(p).((0.0, 20.0)), copy(p);
                      constant_lags = [p[2]])
    solve(prob, MethodOfSteps(Tsit5()); abstol = 1e-14, reltol = 1e-14).u[end]
  end

  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = Calculus.finite_difference_jacobian(test, p)
  fordiff = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.jacobian(test, p)

  @test_broken findiff ≈ fordiff
end

@testset "Hessian" begin
  function test(p)
    prob = DDEProblem(f, [p[5]], h, eltype(p).((0.0, 20.0)), copy(p);
                      constant_lags = [p[2]])
    sol = solve(prob, MethodOfSteps(Tsit5()); abstol = 1e-14, reltol = 1e-14)
    diff = @. sol[1, :] - 10 * exp(-sol.t)
    dot(diff, diff)
  end

  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = Calculus.finite_difference_hessian(test, p)
  fordiff = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.hessian(test, p)

  @test_broken findiff ≈ fordiff
end
