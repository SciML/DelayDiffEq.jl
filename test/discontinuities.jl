using DelayDiffEq, DiffEqProblemLibrary, Base.Test

dde_int = init(prob_dde_2delays(1.0), MethodOfSteps(BS3()))

@test dde_int.tracked_discontinuities == [Discontinuity(0., 0)]

solve!(dde_int)

@test dde_int.tracked_discontinuities ==
    [Discontinuity(t, order) for (t, order) in ((0., 0), (1/5, 1), (1/3, 1),
                                                (2/5, 2), (8/15, 2), (2/3, 2))]
