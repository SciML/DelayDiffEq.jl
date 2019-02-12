include("common.jl")
@testset "Lazy interpolants" begin
    ## simple problems
    @testset "simple problems" begin
        prob_inplace = prob_dde_1delay
        prob_notinplace = prob_dde_1delay_scalar_notinplace
        ts = 0:0.1:10

        # Vern6
        @testset "Vern6" begin
            sol = solve(prob_inplace, MethodOfSteps(Vern6()))

            @test sol.errors[:l∞] < 8.7e-4
            @test sol.errors[:final] < 7.5e-6
            @test sol.errors[:l2] < 5.3e-4

            sol2 = solve(prob_notinplace, MethodOfSteps(Vern6()))

            # fails due to floating point issues on Haswell CPUs: https://github.com/JuliaDiffEq/DelayDiffEq.jl/issues/97
            # @test sol.t ≈ sol2.t && sol[1, :] ≈ sol2.u
            @test sol(ts, idxs=1) ≈ sol2(ts) atol = 1e-8
        end

        # Vern7
        @testset "Vern7" begin
            sol = solve(prob_inplace, MethodOfSteps(Vern7()))

            @test sol.errors[:l∞] < 4.0e-4
            @test sol.errors[:final] < 3.5e-7
            @test sol.errors[:l2] < 1.9e-4

            sol2 = solve(prob_notinplace, MethodOfSteps(Vern7()))

            @test sol.t ≈ sol2.t && sol[1, :] ≈ sol2.u
        end

        # Vern8
        @testset "Vern8" begin
            sol = solve(prob_inplace, MethodOfSteps(Vern8()))

            @test sol.errors[:l∞] < 2.0e-3
            @test sol.errors[:final] < 1.8e-5
            @test sol.errors[:l2] < 9.0e-4

            sol2 = solve(prob_notinplace, MethodOfSteps(Vern8()))

            @test sol.t ≈ sol2.t && sol[1, :] ≈ sol2.u
        end

        # Vern9
        @testset "Vern9" begin
            sol = solve(prob_inplace, MethodOfSteps(Vern9()))

            @test sol.errors[:l∞] < 1.5e-3
            @test sol.errors[:final] < 3.8e-6
            @test sol.errors[:l2] < 6.2e-4

            sol2 = solve(prob_notinplace, MethodOfSteps(Vern9()))

            # fails due to floating point issues on Haswell CPUs: https://github.com/JuliaDiffEq/DelayDiffEq.jl/issues/97
            # @test sol.t ≈ sol2.t && sol[1, :] ≈ sol2.u
            @test sol(ts, idxs=1) ≈ sol2(ts) atol = 1e-12
        end
    end

    # model of Mackey and Glass
    @testset "Mackey and Glass" begin
        prob = prob_dde_mackey

        # Vern6
        solve(prob, MethodOfSteps(Vern6()))

        # Vern7
        solve(prob, MethodOfSteps(Vern7()))

        # Vern8
        solve(prob, MethodOfSteps(Vern8()))

        # Vern9
        solve(prob, MethodOfSteps(Vern9()))
    end
end
