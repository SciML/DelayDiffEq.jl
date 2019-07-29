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
    sol = solve(prob, MethodOfSteps(Tsit5()))
    diff = @. sol[1, :] - 10 * exp(-sol.t)
    dot(diff, diff)
  end

  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = Calculus.finite_difference(test, p)
  fordiff = ForwardDiff.gradient(test, p)

  @test_broken findiff ≈ fordiff
end

@testset "Jacobian" begin
  function test(p)
    prob = DDEProblem(f, [p[5]], h, eltype(p).((0.0, 20.0)), copy(p);
                      constant_lags = [p[2]])
    solve(prob, MethodOfSteps(Tsit5())).u[end]
  end

  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = Calculus.finite_difference_jacobian(test, p)
  fordiff = ForwardDiff.jacobian(test, p)

  @test_broken findiff ≈ fordiff
end

@testset "Hessian" begin
  function test(p)
    prob = DDEProblem(f, [p[5]], h, eltype(p).((0.0, 20.0)), copy(p);
                      constant_lags = [p[2]])
    sol = solve(prob, MethodOfSteps(Tsit5()))
    diff = @. sol[1, :] - 10 * exp(-sol.t)
    dot(diff, diff)
  end

  p = [1.0, 1.0, 1.0, 0.0, 1.0]
  findiff = Calculus.finite_difference_hessian(test, p)
  fordiff = ForwardDiff.hessian(test, p)

  @test_broken findiff ≈ fordiff
end
