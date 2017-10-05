using DelayDiffEq, DiffEqProblemLibrary, Base.Test

## simple problems
prob_inplace = prob_dde_1delay
prob_notinplace = prob_dde_1delay_scalar_notinplace

# Vern6
sol = solve(prob_inplace, MethodOfSteps(Vern6()))

@test sol.errors[:l∞] < 8.7e-4
@test sol.errors[:final] < 7.5e-6
@test sol.errors[:l2] < 5.3e-4

sol2 = solve(prob_notinplace, MethodOfSteps(Vern6()))

@test sol.t == sol2.t && sol[1, :] == sol2.u

# Vern7
sol = solve(prob_inplace, MethodOfSteps(Vern7()))

@test sol.errors[:l∞] < 4.0e-4
@test sol.errors[:final] < 3.5e-7
@test sol.errors[:l2] < 1.9e-4

sol2 = solve(prob_notinplace, MethodOfSteps(Vern7()))

@test sol.t == sol2.t && sol[1, :] == sol2.u

# Vern8

sol = solve(prob_inplace, MethodOfSteps(Vern8()))

@test sol.errors[:l∞] < 2.0e-3
@test sol.errors[:final] < 1.8e-5
@test sol.errors[:l2] < 8.0e-4

sol2 = solve(prob_notinplace, MethodOfSteps(Vern8()))

@test sol.t == sol2.t && sol[1, :] == sol2.u

# Vern9
sol = solve(prob_inplace, MethodOfSteps(Vern9()))

@test sol.errors[:l∞] < 1.5e-3
@test sol.errors[:final] < 3.8e-6
@test sol.errors[:l2] < 6.2e-4

sol2 = solve(prob_notinplace, MethodOfSteps(Vern9()))

@test sol.t == sol2.t && sol[1, :] == sol2.u

# model of Mackey and Glass
prob = prob_dde_mackey

# Vern6
sol = solve(prob, MethodOfSteps(Vern6()))

# Vern7
sol = solve(prob, MethodOfSteps(Vern7()))

# Vern8
sol = solve(prob, MethodOfSteps(Vern8()))

# Vern9
sol = solve(prob, MethodOfSteps(Vern9()))
