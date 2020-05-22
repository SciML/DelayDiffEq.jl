using DelayDiffEq

import FiniteDiff
import ForwardDiff

using LinearAlgebra, Test

# Hutchinson's equation
function f(du, u, h, p, t)
  du[1] = p[2] * u[1] * (1 - h(p, t - p[1])[1] / p[3])
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

  # with delay length estimation
  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = FiniteDiff.finite_difference_gradient(test, p)
  fordiff = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.gradient(test, p)
  @test_broken findiff ≈ fordiff

  # without delay length estimation
  p = [1.0, 1.0, 0.0, 1.0]
  findiff2 = FiniteDiff.finite_difference_gradient(p -> test(vcat(1, p)), p)
  fordiff2 = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.gradient(p -> test(vcat(1, p)), p)
  @test_broken findiff2 ≈ fordiff2

  @test findiff[2:end] ≈ findiff2
end

@testset "Jacobian" begin
  function test(p)
    prob = DDEProblem(f, [p[5]], h, eltype(p).((0.0, 20.0)), copy(p);
                      constant_lags = [p[2]])
    solve(prob, MethodOfSteps(Tsit5()); abstol = 1e-14, reltol = 1e-14).u[end]
  end

  # with delay length estimation
  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = FiniteDiff.finite_difference_jacobian(test, p)
  fordiff = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.jacobian(test, p)
  @test_broken findiff ≈ fordiff

  # without delay length estimation
  p = [1.0, 1.0, 0.0, 1.0]
  findiff2 = FiniteDiff.finite_difference_jacobian(p -> test(vcat(1, p)), p)
  fordiff2 = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.jacobian(p -> test(vcat(1, p)), p)
  @test_broken findiff2 ≈ fordiff2

  @test findiff[:, 2:end] ≈ findiff2
end

@testset "Hessian" begin
  function test(p)
    prob = DDEProblem(f, [p[5]], h, eltype(p).((0.0, 20.0)), copy(p);
                      constant_lags = [p[2]])
    sol = solve(prob, MethodOfSteps(Tsit5()); abstol = 1e-14, reltol = 1e-14)
    diff = @. sol[1, :] - 10 * exp(-sol.t)
    dot(diff, diff)
  end

  # with delay length estimation
  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = FiniteDiff.finite_difference_hessian(test, p)
  fordiff = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.hessian(test, p)
  @test_broken findiff ≈ fordiff

  # without delay length estimation
  p = [1.0, 1.0, 0.0, 1.0]
  findiff2 = FiniteDiff.finite_difference_hessian(p -> test(vcat(1, p)), p)
  fordiff2 = @test_logs (:warn, r"^dt <= dtmin") ForwardDiff.hessian(p -> test(vcat(1, p)), p)
  @test_broken findiff2 ≈ fordiff2

  @test findiff[2:end, 2:end] ≈ findiff2
end
