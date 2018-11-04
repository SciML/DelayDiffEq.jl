include("common.jl")
@testset "Backward integration" begin
    h(p, t) = 1.0
    u₀ = 1.0
    tspan = (2.0, 0.0)

    f(u, h, p,t) = h(p, t + 1)
    # solution on [0,2]
    function f_analytic(u₀, p, t)
        if t < 1
            return (t^2 - 1)/2
        else
            return t - 1
        end
    end
    dde_f = DDEFunction(f, analytic = f_analytic)

    @testset "Without lags" begin
        sol = solve(DDEProblem(dde_f, u₀, h, tspan), MethodOfSteps(RK4()))
        @test sol.errors[:l∞] < 1e-5
    end

    @testset "Constant lags" begin
        # incorrect lags
        prob1 = DDEProblem(dde_f, u₀, h, tspan; constant_lags = [1.0])
        @test_throws ErrorException solve(prob1, MethodOfSteps(Tsit5()))

        prob2 = DDEProblem(dde_f, u₀, h, tspan; constant_lags = [-1.0])
        dde_int = init(prob2, MethodOfSteps(Tsit5()))
        @test dde_int.opts.d_discontinuities.valtree == [Discontinuity(1.0, 2)]
        sol = solve!(dde_int)
        @test sol.errors[:l∞] < 1e-12
        @test dde_int.tracked_discontinuities == [Discontinuity(2.0, 1),
                                                  Discontinuity(1.0, 2)]
    end

    @testset "Dependent lags" begin
        prob = DDEProblem(dde_f, u₀, h, tspan; dependent_lags = ((u, p, t) -> -1.0,))
        dde_int = init(prob, MethodOfSteps(Tsit5()))
        @test isempty(dde_int.opts.d_discontinuities)
        sol = solve!(dde_int)
        @test sol.errors[:l∞] < 1e-12
        @test dde_int.tracked_discontinuities == [Discontinuity(2.0, 1),
                                                  Discontinuity(1.0, 2)]
    end

    @testset "dt and dtmax" begin
        prob = DDEProblem(dde_f, u₀, h, tspan)

        dde_int = init(prob, MethodOfSteps(RK4()); dt = 0.1, dtmax = 0.5)
        @test dde_int.dt == -0.1 && dde_int.opts.dtmax == -0.5
        sol1 = solve!(dde_int)

        sol2 = solve(prob, MethodOfSteps(RK4()); dt = -0.1, dtmax = -0.5)
        @test sol1.t == sol2.t && sol1.u == sol2.u
    end
end
