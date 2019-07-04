include("common.jl")

const T = OrdinaryDiffEq.ODECompositeSolution

const integrator = init(prob_dde_constant_1delay_ip,
                        MethodOfSteps(AutoTsit5(Rosenbrock23())))

@testset "init" begin
  @test integrator.sol isa T
  @test integrator.integrator.sol isa T
end

@testset "solve" begin
  sol1 = solve!(integrator)
  @test sol1 isa T

  # compare integration grid
  sol2 = solve(prob_dde_constant_1delay_ip, MethodOfSteps(Tsit5()))
  @test sol1.t == sol2.t && sol1.u == sol2.u

  # compare interpolation
  @test sol1(0:0.1:1) == sol2(0:0.1:1)
end
