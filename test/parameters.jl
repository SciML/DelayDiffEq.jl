include("common.jl")
@testset "Parameterized functions" begin
    # Test parameterized delayed logistic equation

    # delayed logistic equation
    f_inplace(du, u, h, p, t) = (du[1] = p[1] * u[1] * (1 - h(p, t-1; idxs = 1)))
    f_scalar(u, h, p, t) = p[1] * u * (1 - h(p, t-1))

    # simple history function (simplicity here comes at a price, see below)
    h(p, t; idxs=nothing) = 0.1

    @testset for inplace in (true, false)
        # define problem
        prob = DDEProblem(inplace ? f_inplace : f_scalar,
                          inplace ? [0.1] : 0.1,
                          h, (0.0, 50.0), [0.3]; constant_lags = [1])

        # solve problem with initial parameter:
        # since h(p, 0) = u0 holds only in the scalar case we have to specify
        # initial_order=1 to indicate that the discontinuity at t = 0 is of first
        # order and to obtain the same results in both cases
        sol1 = solve(prob, MethodOfSteps(Tsit5()); initial_order = 1)
        @test length(sol1) == 26
        @test first(sol1(12)) ≈ 0.884 atol=1e-4
        @test first(sol1[end]) ≈ 1 atol=1e-5

        # solve problem with updated parameter
        prob.p[1] = 1.4
        sol2 = solve(prob, MethodOfSteps(Tsit5()); initial_order = 1)
        @test length(sol2) == 52
        @test first(sol2(12)) ≈ 1.127 atol=4e-4
        @test first(sol2[end]) ≈ 0.995 atol=3e-4
    end
end
