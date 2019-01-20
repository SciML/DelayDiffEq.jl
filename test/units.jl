include("common.jl")
using Unitful

@testset "Units" begin
    prob_notinplace =
        DiffEqProblemLibrary.DDEProblemLibrary.build_prob_dde_1delay_long_scalar_notinplace(1.0u"N", 1.0u"s")
    prob_inplace = DiffEqProblemLibrary.DDEProblemLibrary.build_prob_dde_1delay_long(1.0u"N", 1.0u"s")

    @testset for prob in (prob_notinplace, prob_inplace)
        @testset "correct" begin
            # default tolerances
            alg1 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100)
            sol1 = solve(prob, alg1)

            # with correct units
            alg2 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                                 fixedpoint_abstol=1e-6u"N", fixedpoint_reltol=1e-3u"N")
            sol2 = solve(prob, alg2)

            @test sol1.t == sol2.t && sol1.u == sol2.u

            # with correct units as vectors
            if typeof(prob.u0) <: AbstractArray
                alg3 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                                     fixedpoint_abstol=[1e-6u"N"], fixedpoint_reltol=[1e-3u"N"])
                sol3 = solve(prob, alg3)

                @test sol1.t == sol3.t && sol1.u == sol3.u
            end
        end

        @testset "incorrect" begin
            # without units
            alg1 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                                 fixedpoint_abstol=1e-6, fixedpoint_reltol=1e-3)
            @test_throws Unitful.DimensionError solve(prob, alg1)

            # with incorrect units for absolute tolerance
            alg2 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                                fixedpoint_abstol=1e-6u"s", fixedpoint_reltol=1e-3u"N")
            @test_throws Unitful.DimensionError solve(prob, alg2)

            # with incorrect units for relative tolerance
            alg3 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                                fixedpoint_abstol=1e-6u"N", fixedpoint_reltol=1e-3u"s")
            @test_throws Unitful.DimensionError solve(prob, alg3)
        end
    end
end
