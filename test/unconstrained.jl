include("common.jl")
@testset "Unconstrained time step" begin
    # Check that numerical solutions approximate analytical solutions,
    # independent of problem structure
    @testset "standard history" begin
        # standard algorithm
        alg = MethodOfSteps(BS3(); constrained=false)

        # different tolerances
        alg1 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                             fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)
        alg2 = MethodOfSteps(DP8(), constrained=false, max_fixedpoint_iters=10,
                             fixedpoint_abstol=1e-8, fixedpoint_reltol=1e-10)
        alg3 = MethodOfSteps(Tsit5(), constrained=true, max_fixedpoint_iters=10,
                             fixedpoint_abstol=1e-8, fixedpoint_reltol=1e-10)
        alg4 = MethodOfSteps(DP5(), constrained=false, max_fixedpoint_iters=100,
                             fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)

        ## Single constant delay
        @testset "single constant delay" begin
            @testset "short time span" begin
                ### Not in-place function with scalar history function
                prob = prob_dde_1delay_scalar_notinplace
                sol = solve(prob, alg)

                @test sol.errors[:l∞] < 3.0e-5
                @test sol.errors[:final] < 2.1e-5
                @test sol.errors[:l2] < 1.2e-5

                ### Not in-place function with vectorized history function
                prob = prob_dde_1delay_notinplace
                sol2 = solve(prob, alg)

                @test sol.t ≈ sol2.t && sol.u ≈ sol2[1, :]

                ### In-place function
                prob = prob_dde_1delay
                sol2 = solve(prob, alg)

                @test sol.t ≈ sol2.t && sol.u ≈ sol2[1, :]
            end

            @testset "long time span" begin
                prob = prob_dde_1delay_long_scalar_notinplace

                sol1 = solve(prob, alg1)
                sol2 = solve(prob, alg2)
                sol3 = solve(prob, alg3)
                sol4 = solve(prob, alg4)

                @test abs(sol1[end] - sol2[end]) < 5.3e-4
                @test abs(sol1[end] - sol3[end]) < 1.1e-8
                @test abs(sol1[end] - sol4[end]) < 2.9e-4
            end
        end

        ## Two constant delays
        @testset "two constant delays" begin
            @testset "short time span" begin
                ### Not in-place function with scalar history function
                prob = prob_dde_2delays_scalar_notinplace
                sol = solve(prob, alg)

                @test sol.errors[:l∞] < 1.9e-6
                @test sol.errors[:final] < 1.2e-6
                @test sol.errors[:l2] < 1.0e-6

                ### Not in-place function with vectorized history function
                prob = prob_dde_2delays_notinplace
                sol2 = solve(prob, alg)

                @test sol.t ≈ sol2.t && sol.u ≈ sol2[1, :]

                ### In-place function
                prob = prob_dde_2delays
                sol2 = solve(prob, alg)

                @test sol.t ≈ sol2.t && sol.u ≈ sol2[1, :]
            end

            @testset "long time span" begin
                prob = prob_dde_2delays_long_scalar_notinplace

                sol1 = solve(prob, alg1)
                # sol2 = solve(prob, alg) # aborted because of dt <= dtmin
                sol3 = solve(prob, alg3)
                sol4 = solve(prob, alg4)

                # relaxed tests to prevent floating point issues
                @test abs(sol1[end] - sol3[end]) < 7.9e-13 # 7.9e-15
                @test abs(sol1[end] - sol4[end]) < 1.6e-12 # 1.6e-14
            end
        end
    end

    ## Non-standard history functions
    @testset "non-standard history" begin
        alg = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                            fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)

        @testset "idxs" begin
            function f(du,u,h,p,t)
                du[1] = -h(p, t-0.2;idxs=1) + u[1]
            end
            h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : [0.0]

            prob = DDEProblem(f, [1.0], h, (0.0, 100.), constant_lags = [0.2])
            solve(prob, alg)
        end

        @testset "in-place" begin
            function f(du,u,h,p,t)
                h(du, p, t-0.2)
                du[1] = -du[1] + u[1]
            end
            h(val, p, t) = (val .= 0.0)

            prob = DDEProblem(f, [1.0], h, (0.0, 100.0), constant_lags = [0.2])
            solve(prob, alg)
        end
    end
end
