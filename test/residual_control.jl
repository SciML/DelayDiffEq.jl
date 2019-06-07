include("common.jl")
@testset "Residual control" begin
    alg = MethodOfSteps(RK4(); constrained=false)

    # reference solution with delays specified
    @testset "reference" begin
        prob = prob_dde_1delay_scalar_notinplace
        sol = solve(prob, alg)

        @test sol.errors[:l∞] < 5.6e-5
        @test sol.errors[:final] < 1.8e-6
        @test sol.errors[:l2] < 2.0e-5
    end

    # problem without delays specified
    prob = DDEProblem(prob.f,prob.u0,prob.h,prob.tspan)

    # solutions with residual control
    @testset "residual control" begin
        sol = solve(prob, alg)

        @test sol.errors[:l∞] < 1.8e-4
        @test sol.errors[:final] < 4.1e-6
        @test sol.errors[:l2] < 9.0e-5

        sol = solve(prob, alg, abstol=1e-9,reltol=1e-6)

        @test sol.errors[:l∞] < 1.5e-7
        @test sol.errors[:final] < 4.1e-9
        @test sol.errors[:l2] < 7.5e-8

        sol = solve(prob, alg, abstol=1e-13,reltol=1e-13)

        @test sol.errors[:l∞] < 7.0e-11
        @test sol.errors[:final] < 1.1e-11
        @test sol.errors[:l2] < 9.3e-12
    end

    ######## Now show that non-residual control is worse
    # solutions without residual control
    @testset "non-residual control" begin
        alg = MethodOfSteps(OwrenZen5(); constrained=false)
        sol = solve(prob, alg)

        @test sol.errors[:l∞] > 1e-1
        @test sol.errors[:final] > 1e-3
        @test sol.errors[:l2] > 4e-2

        alg = MethodOfSteps(OwrenZen5(); constrained=true)
        sol = solve(prob, alg)

        @test sol.errors[:l∞] > 1e-1
        @test sol.errors[:final] > 1e-3
        @test sol.errors[:l2] > 4e-2

        sol = solve(prob, alg,abstol=1e-13,reltol=1e-13)

        @test sol.errors[:l∞] > 1e-1
        @test sol.errors[:final] > 1e-3
        @test sol.errors[:l2] > 5e-2
    end
end
