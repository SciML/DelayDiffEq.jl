@testset "solutions" begin
    (; default, anderson) = compare_anderson()

    # Check success and desired properties
    for sol in (default, anderson)
        @test sol.retcode === :Success
        @test sol.prob ===
              prob_dde_constant_2delays_long_scalar
    end

    # Compute reference solution with low tolerances
    alg = MethodOfSteps(Vern9();
                        fpsolve = NLFunctional(;
                                               max_iter = 1000))
    testsol = TestSolution(solve(default.prob, alg;
                                 reltol = 1e-14,
                                 abstol = 1e-14))

    # Check errors of solutions with and without Anderson acceleration
    # Note that the range of the errors is very similar
    default2 = appxtrue(default, testsol)
    @test default2.errors[:L2] < 3e-5
    @test default2.errors[:final] < 4e-9
    @test default2.errors[:L∞] < 3e-4

    anderson2 = appxtrue(anderson, testsol)
    @test anderson2.errors[:L2] < 3e-5
    @test anderson2.errors[:final] < 6e-10
    @test anderson2.errors[:L∞] < 3e-4
end

@testset "statistics" begin
    (; default, anderson) = compare_anderson()

    # Statistics with default algorithm
    defaultstats = default.destats
    @test defaultstats.nf == 5517
    @test defaultstats.nfpiter == 720
    @test defaultstats.nfpconvfail == 58

    # Statistics with Anderson acceleration
    andersonstats = anderson.destats
    @test andersonstats.nf == 2829
    @test andersonstats.nfpiter == 308
    @test andersonstats.nfpconvfail == 30
end
