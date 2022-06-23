using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using Test

DDEProblemLibrary.importddeproblems()

# simple problems
@testset "simple problems" begin
    prob_ip = DDEProblemLibrary.prob_dde_constant_1delay_ip
    prob_scalar = DDEProblemLibrary.prob_dde_constant_1delay_scalar
    ts = 0:0.1:10

    # Vern6
    println("Vern6")
    @testset "Vern6" begin
        alg = MethodOfSteps(Vern6())
        sol_ip = solve(prob_ip, alg)

        @test sol_ip.errors[:l∞] < 8.7e-4
        @test sol_ip.errors[:final] < 3.9e-6
        @test sol_ip.errors[:l2] < 5.4e-4

        sol_scalar = solve(prob_scalar, alg)

        # fails due to floating point issues
        if Sys.WORD_SIZE == 32
            @test_broken sol_ip(ts, idxs = 1) ≈ sol_scalar(ts)
        else
            @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts)
        end
    end

    # Vern7
    println("Vern7")
    @testset "Vern7" begin
        alg = MethodOfSteps(Vern7())
        sol_ip = solve(prob_ip, alg)

        @test sol_ip.errors[:l∞] < 4.0e-4
        @test sol_ip.errors[:final] < 3.5e-7
        @test sol_ip.errors[:l2] < 1.9e-4

        sol_scalar = solve(prob_scalar, alg)

        @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts)
    end

    # Vern8
    println("Vern8")
    @testset "Vern8" begin
        alg = MethodOfSteps(Vern8())
        sol_ip = solve(prob_ip, alg)

        @test sol_ip.errors[:l∞] < 2.0e-3
        @test sol_ip.errors[:final] < 1.8e-5
        @test sol_ip.errors[:l2] < 8.8e-4

        sol_scalar = solve(prob_scalar, alg)

        @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts)
    end

    # Vern9
    println("Vern9")
    @testset "Vern9" begin
        alg = MethodOfSteps(Vern9())
        sol_ip = solve(prob_ip, alg)

        @test sol_ip.errors[:l∞] < 1.5e-3
        @test sol_ip.errors[:final] < 3.8e-6
        @test sol_ip.errors[:l2] < 6.5e-4

        sol_scalar = solve(prob_scalar, alg)

        @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts)
    end
end

# model of Mackey and Glass
println("Mackey and Glass")
@testset "Mackey and Glass" begin
    prob = DDEProblemLibrary.prob_dde_DDETST_A1

    # Vern6
    solve(prob, MethodOfSteps(Vern6()))

    # Vern7
    solve(prob, MethodOfSteps(Vern7()))

    # Vern8
    solve(prob, MethodOfSteps(Vern8()))

    # Vern9
    solve(prob, MethodOfSteps(Vern9()))
end
