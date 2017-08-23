using DelayDiffEq, DiffEqProblemLibrary, Base.Test

alg = MethodOfSteps(RK4(); constrained=false)
u₀ = 1.0

prob = prob_dde_1delay_scalar_notinplace(u₀)
sol = solve(prob, alg)

@test sol.errors[:l∞] < 1e-4
@test sol.errors[:final] < 1e-4
@test sol.errors[:l2] < 1e-4

prob2 = deepcopy(prob)
typeof(prob2) <: DDEConstantLagProblem ? prob2.lags=Rational{Int}[] :
                 prob2.constant_lags=Rational{Int}[]
sol2 = solve(prob2, alg)

@test sol2.errors[:l∞] < 1e-1
@test sol2.errors[:final] < 1e-1
@test sol2.errors[:l2] < 1e-1

sol3 = solve(prob2, alg, abstol=1e-9,reltol=1e-6)

@test sol3.errors[:l∞] < 1e-1
@test sol3.errors[:final] < 1e-3
@test sol3.errors[:l2] < 1e-1

sol4 = solve(prob2, alg, abstol=1e-13,reltol=1e-13)

@test sol4.errors[:l∞] < 1e-1
@test sol4.errors[:final] < 1e-3
@test sol4.errors[:l2] < 1e-1

######## Now show that non-residual control is worse

alg = MethodOfSteps(OwrenZen5(); constrained=false)

sol4 = solve(prob2, alg)

@test sol4.errors[:l∞] > 1e-1
@test sol4.errors[:final] > 1e-3
@test sol4.errors[:l2] > 1e-2

alg = MethodOfSteps(OwrenZen5(); constrained=true)

sol5 = solve(prob2, alg)

@test sol5.errors[:l∞] > 1e-1
@test sol5.errors[:final] > 1e-3
@test sol5.errors[:l2] > 1e-2

sol5 = solve(prob2, alg,abstol=1e-13,reltol=1e-13)

@test sol5.errors[:l∞] > 1e-1
@test sol5.errors[:final] > 1e-3
@test sol5.errors[:l2] > 1e-2
