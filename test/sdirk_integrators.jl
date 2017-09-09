using DelayDiffEq, DiffEqProblemLibrary, Base.Test

prob = prob_dde_1delay

# ODE algorithms
algs = [GenericImplicitEuler(), GenericTrapezoid(),
        ImplicitEuler(), ImplicitMidpoint(), Trapezoid(),
        TRBDF2(), SDIRK2(), SSPSDIRK2(),
        Kvaerno3(), KenCarp3(),
        Cash4(), Hairer4(), Hairer42(), Kvaerno4(), KenCarp4(),
        Kvaerno5(), KenCarp5()]

names = ["GenericImplicitEuler", "GenericTrapezoid",
         "ImplicitEuler", "ImplicitMidpoint", "Trapezoid",
         "TRBDF2", "SDIRK2", "SSPSDIRK2",
         "Kvaerno3", "KenCarp3",
         "Cash4", "Hairer4", "Hairer42", "Kvaerno4", "KenCarp4",
         "Kvaerno5", "KenCarp5"]

for (alg, name) in zip(algs, names)
    print("testing ", name, "... ")
    step_alg = MethodOfSteps(alg)
    solve(prob, step_alg)
    @time solve(prob, step_alg)
end
