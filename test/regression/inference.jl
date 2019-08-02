using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using Test

DDEProblemLibrary.importddeproblems()

@testset "init" begin
    prob = DDEProblemLibrary.prob_dde_constant_1delay_ip
    prob_scalar = DDEProblemLibrary.prob_dde_constant_1delay_scalar

    inferred = [BS3(), Tsit5(), RK4(), Vern6()]
    for alg in inferred
        ddealg = MethodOfSteps(alg)

        @inferred init(prob, ddealg)
        @inferred init(prob_scalar, ddealg)
    end

    notinferred = [SDIRK2(), TRBDF2(), KenCarp4(), Rosenbrock23(), Rodas4()]
    for alg in notinferred
        ddealg = MethodOfSteps(alg)

        @test_broken @inferred init(prob, ddealg)
        @test_broken @inferred init(prob_scalar, ddealg)
    end
end