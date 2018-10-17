include("common.jl")
@testset "Constrained time step" begin
    # Check that numerical solutions approximate analytical solutions,
    # independent of problem structure

    alg = MethodOfSteps(BS3(); constrained=true)

    # Single constant delay
    @testset "single constant delay" begin
        ## Not in-place function with scalar history function
        prob = prob_dde_1delay_scalar_notinplace
        dde_int = init(prob, alg; dt=0.1)
        sol = solve!(dde_int)

        @test sol.errors[:l∞] < 3e-5
        @test sol.errors[:final] < 2.1e-5
        @test sol.errors[:l2] < 1.3e-5

        ## Not in-place function with vectorized history function
        prob = prob_dde_1delay_notinplace
        dde_int = init(prob, alg; dt=0.1)
        sol2 = solve!(dde_int)

        @test sol.t ≈ sol2.t && sol.u ≈ sol2[1, :]

        ## In-place function
        prob = prob_dde_1delay
        dde_int = init(prob, alg; dt=0.1)
        sol2 = solve!(dde_int)

        @test sol.t == sol2.t && sol.u == sol2[1, :]
    end

    # Two constant delays
    @testset "two constant delays" begin
        ## Not in-place function with scalar history function
        prob = prob_dde_2delays_scalar_notinplace
        dde_int = init(prob, alg; dt=0.1)
        sol = solve!(dde_int)

        @test sol.errors[:l∞] < 4.1e-6
        @test sol.errors[:final] < 1.5e-6
        @test sol.errors[:l2] < 2.3e-6

        ## Not in-place function with vectorized history function
        prob = prob_dde_2delays_notinplace
        dde_int = init(prob, alg; dt=0.1)
        sol2 = solve!(dde_int)

        @test sol.t ≈ sol2.t && sol.u ≈ sol2[1, :]

        ## In-place function
        prob = prob_dde_2delays
        dde_int = init(prob, alg; dt=0.1)
        sol2 = solve!(dde_int)

        @test sol.t == sol2.t && sol.u == sol2[1, :]
    end
end
