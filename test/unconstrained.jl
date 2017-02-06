using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, Base.Test

lags = [.2]

f = function (t,u,h,du)
  du[1] = -h(t-.2)[1] + u[1]
end
h = (t) -> [0.0]

f = function (t,u,h)
  out = -h(t-.2) + u
end
h = (t) -> 0.0

prob = ConstantLagDDEProblem(f,h,1.0,lags,(0.0,100.0))

alg1 = MethodOfSteps(Tsit5(),constrained=false,max_picard_iters=100,picardabstol=1e-12,picardreltol=1e-12)
sol1 = solve(prob,alg1)

alg2 = MethodOfSteps(DP8(),constrained=false,max_picard_iters=10,picardabstol=1e-8,picardreltol=1e-10)
sol2 = solve(prob,alg2)

alg3 = MethodOfSteps(Tsit5(),constrained=true,max_picard_iters=10,picardabstol=1e-8,picardreltol=1e-10)
sol3 = solve(prob,alg3)

alg4 = MethodOfSteps(DP5(),constrained=false,max_picard_iters=100,picardabstol=1e-12,picardreltol=1e-12)
sol4 = solve(prob,alg4)

@test abs(sol1[end] - sol2[end]) < 1e-3
@test abs(sol1[end] - sol3[end]) < 1e-3
@test abs(sol1[end] - sol4[end]) < 1e-3

lags = [1//3,1//5]
f = function (t,u,h)
  - h(t-1/3) - h(t-1/5)
end
h = (t) -> 0.0

prob = ConstantLagDDEProblem(f,h,1.0,lags,(0.0,10.0);iip=DiffEqBase.isinplace(f,4))

alg1 = MethodOfSteps(Tsit5(),constrained=false,max_picard_iters=100,picardabstol=1e-12,picardreltol=1e-12)
sol1 = solve(prob,alg1)

alg2 = MethodOfSteps(DP8(),constrained=false,max_picard_iters=10,picardabstol=1e-8,picardreltol=1e-10)
sol2 = solve(prob,alg2)

alg3 = MethodOfSteps(Tsit5(),constrained=true,max_picard_iters=10,picardabstol=1e-8,picardreltol=1e-10)
sol3 = solve(prob,alg3)

alg4 = MethodOfSteps(DP5(),constrained=false,max_picard_iters=100,picardabstol=1e-12,picardreltol=1e-12)
sol4 = solve(prob,alg4)

@test abs(sol1[end] - sol2[end]) < 1e-3
@test abs(sol1[end] - sol3[end]) < 1e-3
@test abs(sol1[end] - sol4[end]) < 1e-3

println("Standard tests complete. Onto idxs tests")

## Idxs

f = function (t,u,h,du)
  du[1] = -h(t-.2,Val{0},1) + u[1]
end

h = function (t,idxs=nothing)
  if typeof(idxs) <: Void
    return [0.0]
  else
    return 0.0
  end
end

prob = ConstantLagDDEProblem(f,h,[1.0],lags,(0.0,100.0))

alg1 = MethodOfSteps(Tsit5(),constrained=false,max_picard_iters=100,picardabstol=1e-12,picardreltol=1e-12)
@time sol1 = solve(prob,alg1)

f = function (t,u,h,du)
  h(du,t-.2)
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


prob = ConstantLagDDEProblem(f,h,[1.0],lags,(0.0,100.0))

alg1 = MethodOfSteps(Tsit5(),constrained=false,max_picard_iters=100,picardabstol=1e-12,picardreltol=1e-12)
@time sol1 = solve(prob,alg1)
