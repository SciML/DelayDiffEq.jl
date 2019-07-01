include("common.jl")
using DiffEqDevTools, DiffEqCallbacks

const prob = prob_dde_constant_1delay_scalar
const alg = MethodOfSteps(Tsit5(); constrained=false)

# continuous callback
@testset "continuous" begin
  cb = ContinuousCallback(
    (u, t, integrator) -> t - 2.60, # Event when event_f(t,u,k) == 0
    integrator -> (integrator.u = - integrator.u))

  sol1 = solve(prob, alg, callback=cb)
  sol2 = solve(prob, alg, callback=cb, dtmax=0.01)
  sol3 = appxtrue(sol1, sol2)

  @test sol3.errors[:L2] < 4.1e-3
  @test sol3.errors[:L∞] < 1.3e-2
end

# discrete callback
@testset "discrete" begin
  cb = AutoAbstol()

  sol1 = solve(prob, alg, callback=cb)
  sol2 = solve(prob, alg, callback=cb, dtmax=0.01)
  sol3 = appxtrue(sol1, sol2)

  @test sol3.errors[:L2] < 1.4e-3
  @test sol3.errors[:L∞] < 4.1e-3
end

@testset "save discontinuity" begin
  f(du, u, h, p, t) = (du .= 0)
  prob = DDEProblem(f, [0.0], nothing, (0.0, 1.0))

  condition(u, t, integrator) = t == 0.5
  global iter
  function affect!(integrator) integrator.u[1] += 100; global iter = integrator.iter end
  cb = DiscreteCallback(condition, affect!)
  sol = solve(prob, MethodOfSteps(Tsit5()), callback=cb, tstops=[0.5])
  @test sol.t[iter+1] == sol.t[iter+2]
  @test sol.u[iter+1] == [0.0]
  @test sol.u[iter+2] != [0.0]
end
