using DelayDiffEq, DiffEqProblemLibrary, Base.Test

alg = MethodOfSteps(BS3())

dde_int = init(prob_dde_1delay, alg)

sol = solve!(dde_int)

@test sol.errors[:l∞] < 3.7e-5
@test sol.errors[:final] < 2.0e-5
@test sol.errors[:l2] < 1.5e-5

# constant delay specified as function
prob2 = DDEProblem(DiffEqProblemLibrary.f_1delay, t -> [0.0], [1.0], (0., 10.), [],
                   [(t, u) -> 1])

dde_int2 = init(prob2, alg)
sol2 = solve!(dde_int2)

@test dde_int.tracked_discontinuities == dde_int2.tracked_discontinuities

@test sol2.errors[:l∞] < 3.1e-5
@test sol2.errors[:final] < 1.9e-5
@test sol2.errors[:l2] < 1.3e-5

# without any delays specified is worse
prob3 = DDEProblem(DiffEqProblemLibrary.f_1delay, t -> [0.0], [1.0], (0., 10.), [])

dde_int3 = init(prob3, alg)
sol3 = solve!(dde_int3)

@test sol3.errors[:l∞] > 1e-3
@test sol3.errors[:final] > 4e-5
@test sol3.errors[:l2] > 7e-4
