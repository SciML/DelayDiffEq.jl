using DelayDiffEq, DiffEqProblemLibrary, Base.Test

alg = MethodOfSteps(BS3())

dde_int = init(prob_dde_1delay, alg)

sol = solve!(dde_int)

@test sol.errors[:l∞] < 3.0e-5
@test sol.errors[:final] < 2.1e-5
@test sol.errors[:l2] < 1.2e-5

# constant delay specified as function
prob2 = DDEProblem(DiffEqProblemLibrary.f_1delay, t -> [0.0], [1.0], (0., 10.), nothing, [],
                   [(u,p,t) -> 1])

dde_int2 = init(prob2, alg)
sol2 = solve!(dde_int2)

@test dde_int.tracked_discontinuities == dde_int2.tracked_discontinuities

# with nothing
prob2_nothing = DDEProblem(DiffEqProblemLibrary.f_1delay, t -> [0.0], [1.0], (0., 10.), nothing, nothing,
                   [(u,p,t) -> 1])

dde_int2_nothing = init(prob2_nothing, alg)
sol2_nothing = solve!(dde_int2_nothing)

@test dde_int.tracked_discontinuities == dde_int2_nothing.tracked_discontinuities
@test sol2.u == sol2_nothing.u && sol2.t == sol2_nothing.t

# worse than results above with constant delays specified as scalars
@test sol2.errors[:l∞] < 4.2e-5
@test sol2.errors[:final] < 2.2e-5
@test sol2.errors[:l2] < 1.7e-5

# simple convergence tests
sol3 = solve(prob2, alg, abstol=1e-9, reltol=1e-6)

@test sol3.errors[:l∞] < 7.5e-8
@test sol3.errors[:final] < 4.6e-8
@test sol3.errors[:l2] < 3.9e-8

sol4 = solve(prob2, alg, abstol=1e-13, reltol=1e-13)

@test sol4.errors[:l∞] < 6.9e-11
@test sol4.errors[:final] < 1.1e-11
@test sol4.errors[:l2] < 6.8e-12 # 6.7e-12, relaxed for Win32

# without any delays specified is worse
prob3 = DDEProblem(DiffEqProblemLibrary.f_1delay, t -> [0.0], [1.0], (0., 10.), nothing, [])

dde_int3 = init(prob3, alg)
sol3 = solve!(dde_int3)

@test sol3.errors[:l∞] > 1e-3
@test sol3.errors[:final] > 4e-5
@test sol3.errors[:l2] > 7e-4
