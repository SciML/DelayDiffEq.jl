include("common.jl")
@testset "Unique times" begin
    prob = prob_dde_1delay_long

    @testset for constrained in (false, true)
        alg = MethodOfSteps(Tsit5(), constrained=constrained)
        sol = solve(prob, alg)

        @test allunique(sol.t)
    end
end
