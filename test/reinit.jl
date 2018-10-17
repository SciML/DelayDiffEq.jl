include("common.jl")
using RecursiveArrayTools

@testset "Reinitialization" begin
    alg = MethodOfSteps(BS3(); constrained=false)

    @testset for inplace in (true, false)
        prob = inplace ? prob_dde_1delay : prob_dde_1delay_scalar_notinplace

        @testset "integrator" begin
            integrator = init(prob, alg, dt= 0.01)
            solve!(integrator)

            u = recursivecopy(integrator.sol.u)
            t = copy(integrator.sol.t)

            reinit!(integrator)
            integrator.dt = 0.01
            solve!(integrator)

            @test u == integrator.sol.u
            @test t == integrator.sol.t
        end

        @testset "solution" begin
            integrator = init(prob, alg, dt= 0.01, tstops = [0.5], saveat = [0.33])
            sol = solve!(integrator)

            u = recursivecopy(sol.u)
            t = copy(sol.t)

            reinit!(integrator)
            integrator.dt = 0.01
            sol = solve!(integrator)

            @test u == sol.u
            @test t == sol.t
        end
    end
end
