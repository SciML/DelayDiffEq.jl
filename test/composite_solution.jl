include("common.jl")
@testset "Composite solutions" begin
    T = OrdinaryDiffEq.ODECompositeSolution

    integrator = init(prob_dde_1delay,
                      MethodOfSteps(AutoTsit5(Rosenbrock23())))

    @testset "init" begin
        @test typeof(integrator.sol) <: T
        @test typeof(integrator.integrator.sol) <: T
    end

    @testset "solve" begin
        sol1 = solve!(integrator)
        @test typeof(sol1) <: T

        # compare integration grid
        sol2 = solve(prob_dde_1delay, MethodOfSteps(Tsit5()))
        @test sol1.t == sol2.t && sol1.u == sol2.u

        # compare interpolation
        @test sol1(0:0.1:1) == sol2(0:0.1:1)
    end
end
