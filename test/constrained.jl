using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

# Check that numerical solutions approximate analytical solutions,
# independent of problem structure

alg = MethodOfSteps(BS3(); constrained=true)
u₀ = 1.0

# Single constant delay

## Not in-place function with scalar history function

prob = prob_dde_1delay_scalar_notinplace(u₀)
dde_int = init(prob, alg; dt=0.1)
sol = solve!(dde_int)

@test sol.errors[:l∞] < 3e-5
@test sol.errors[:final] < 2e-5
@test sol.errors[:l2] < 2e-5

## Not in-place function with vectorized history function

prob = prob_dde_1delay_notinplace(u₀)
dde_int = init(prob, alg; dt=0.1)
sol2 = solve!(dde_int)

@test sol.t == sol2.t && sol.u == sol2[1, :]

## In-place function

prob = prob_dde_1delay(u₀)
dde_int = init(prob, alg; dt=0.1)
sol2 = solve!(dde_int)

@test sol.t == sol2.t && sol.u == sol2[1, :]

# Two constant delays

## Not in-place function with scalar history function

prob = prob_dde_2delays_scalar_notinplace(u₀)
dde_int = init(prob, alg; dt=0.1)
sol = solve!(dde_int)

@test sol.errors[:l∞] < 5e-6
@test sol.errors[:final] < 2e-6
@test sol.errors[:l2] < 3e-6

## Not in-place function with vectorized history function

prob = prob_dde_2delays_notinplace(u₀)
dde_int = init(prob, alg; dt=0.1)
sol2 = solve!(dde_int)

@test sol.t == sol2.t && sol.u == sol2[1, :]

## In-place function

prob = prob_dde_2delays(u₀)
dde_int = init(prob, alg; dt=0.1)
sol2 = solve!(dde_int)

@test sol.t == sol2.t && sol.u == sol2[1, :]
