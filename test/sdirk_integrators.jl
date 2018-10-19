include("common.jl")
@testset "SDIRK integrators" begin
    prob_inplace = prob_dde_1delay
    prob_notinplace = prob_dde_1delay_scalar_notinplace

    # ODE algorithms
    algs = [GenericImplicitEuler(), GenericTrapezoid(),
            ImplicitEuler(), ImplicitMidpoint(), Trapezoid(),
            TRBDF2(), SDIRK2(), SSPSDIRK2(),
            Kvaerno3(), KenCarp3(),
            Cash4(), Hairer4(), Hairer42(), Kvaerno4(), KenCarp4(),
            Kvaerno5(), KenCarp5()]

    @testset for alg in algs
        @show alg
        stepsalg = MethodOfSteps(alg)
        solve(prob_inplace, stepsalg)
        solve(prob_notinplace, stepsalg)
    end
end
