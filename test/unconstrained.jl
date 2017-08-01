using DelayDiffEq, OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

# Check that numerical solutions approximate analytical solutions,
# independent of problem structure

alg = MethodOfSteps(BS3(); constrained=false)
u₀ = 1.0

## Single constant delay

### Not in-place function with scalar history function

prob = prob_dde_1delay_scalar_notinplace(u₀)
sol = solve(prob, alg)

@test sol.errors[:l∞] < 3.7e-5
@test sol.errors[:final] < 2.0e-5
@test sol.errors[:l2] < 1.5e-5

### Not in-place function with vectorized history function

prob = prob_dde_1delay_notinplace(u₀)
sol2 = solve(prob, alg)

@test sol.t == sol2.t && sol.u == sol2[1, :]

### In-place function

prob = prob_dde_1delay(u₀)
sol2 = solve(prob, alg)

@test sol.t == sol2.t && sol.u == sol2[1, :]

## Two constant delays

### Not in-place function with scalar history function

prob = prob_dde_2delays_scalar_notinplace(u₀)
sol = solve(prob, alg)

@test sol.errors[:l∞] < 1.9e-6
@test sol.errors[:final] < 1.2e-6
@test sol.errors[:l2] < 1.1e-6

### Not in-place function with vectorized history function

prob = prob_dde_2delays_notinplace(u₀)
sol2 = solve(prob, alg)

@test sol.t == sol2.t && sol.u == sol2[1, :]

### In-place function

prob = prob_dde_2delays(u₀)
sol2 = solve(prob, alg)

@test sol.t == sol2.t && sol.u == sol2[1, :]

# Problems with long time span

## Single constant delay

prob = prob_dde_1delay_long_scalar_notinplace

alg1 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                     fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)
sol1 = solve(prob, alg1)

alg2 = MethodOfSteps(DP8(), constrained=false, max_fixedpoint_iters=10,
                     fixedpoint_abstol=1e-8, fixedpoint_reltol=1e-10)
sol2 = solve(prob, alg2)

alg3 = MethodOfSteps(Tsit5(), constrained=true, max_fixedpoint_iters=10,
                     fixedpoint_abstol=1e-8, fixedpoint_reltol=1e-10)
sol3 = solve(prob, alg3)

alg4 = MethodOfSteps(DP5(), constrained=false, max_fixedpoint_iters=100,
                     fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)
sol4 = solve(prob, alg4)

@test abs(sol1[end] - sol2[end]) < 5.3e-4
@test abs(sol1[end] - sol3[end]) < 4.2e-8
@test abs(sol1[end] - sol4[end]) < 2.9e-4

## Two constant delays

prob = prob_dde_2delays_long_scalar_notinplace

alg1 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                     fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)
sol1 = solve(prob, alg1)

alg2 = MethodOfSteps(DP8(), constrained=false, max_fixedpoint_iters=10,
                     fixedpoint_abstol=1e-8, fixedpoint_reltol=1e-10)
sol2 = solve(prob, alg2)

alg3 = MethodOfSteps(Tsit5(), constrained=true, max_fixedpoint_iters=10,
                     fixedpoint_abstol=1e-8, fixedpoint_reltol=1e-10)
sol3 = solve(prob, alg3)

alg4 = MethodOfSteps(DP5(), constrained=false, max_fixedpoint_iters=100,
                     fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)
sol4 = solve(prob, alg4)

@test abs(sol1[end] - sol2[end]) < 2.5e-11
@test abs(sol1[end] - sol3[end]) < 1.3e-14
@test abs(sol1[end] - sol4[end]) < 7.6e-15

println("Standard tests complete. Onto idxs tests")

# Idxs

f = function (t,u,h,du)
  du[1] = -h(t-0.2, Val{0}, 1) + u[1]
end

h = function (t,idxs=nothing)
  if typeof(idxs) <: Void
    return [0.0]
  else
    return 0.0
  end
end

prob = ConstantLagDDEProblem(f, h, [1.0], [0.2], (0.0, 100.0))

alg1 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                     fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)
@time sol1 = solve(prob, alg1)

f = function (t,u,h,du)
  h(du, t-0.2)
  du[1] = -du[1]
  du[1] += u[1]
end
h = function (t,idxs=nothing)
  if typeof(idxs) <: Void
    return [0.0]
  else
    return 0.0
  end
end
h = function (out,t,idxs=nothing)
  out[1] = 0.0
end

prob = ConstantLagDDEProblem(f, h, [1.0], [0.2], (0.0, 100.0))

alg1 = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                     fixedpoint_abstol=1e-12, fixedpoint_reltol=1e-12)
@time sol1 = solve(prob, alg1)
