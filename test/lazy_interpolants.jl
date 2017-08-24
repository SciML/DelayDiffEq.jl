using DelayDiffEq, DiffEqProblemLibrary, Base.Test

# in-place problem
prob = prob_dde_1delay(1.0)

# Vern7
sol = solve(prob, MethodOfSteps(Vern7()))

@test sol.errors[:l∞] < 4.0e-4
@test sol.errors[:final] < 3.5e-7
@test sol.errors[:l2] < 1.8e-4

# Vern9
sol2 = solve(prob, MethodOfSteps(Vern9()))

@test sol2.errors[:l∞] < 1.5e-3
@test sol2.errors[:final] < 3.8e-6
@test sol2.errors[:l2] < 6.2e-4

# not in-place problem
prob = prob_dde_1delay_scalar_notinplace(1.0)

# Vern7
sol_scalar = solve(prob, MethodOfSteps(Vern7()))

@test sol.t == sol_scalar.t && sol[1, :] == sol_scalar.u

# Vern9
sol2_scalar = solve(prob, MethodOfSteps(Vern9()))

@test sol2.t == sol2_scalar.t && sol2[1, :] == sol2_scalar.u
