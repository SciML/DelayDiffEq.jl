using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, Base.Test

const τ = .2
lags = [τ]

f = function (t,u,h,du)
  du[1] = -h(t-τ)[1] + u[1]
end
h = (t) -> [0.0]

f = function (t,u,h)
  out = -h(t-τ) + u
end
h = (t) -> 0.0


prob = DDEProblem(f,h,1.0,lags,(0.0,1000.0);iip=DiffEqBase.isinplace(f,4))

alg1 = MethodOfSteps(Tsit5(),constrained=false,max_picard_iters=100,picardabstol=1e-12,picardreltol=1e-12)
sol1 = solve(prob,alg1)
using Plots; plot(sol1)

alg2 = MethodOfSteps(DP8(),constrained=false,max_picard_iters=10,picardabstol=1e-8,picardreltol=1e-10)
sol2 = solve(prob,alg2)
using Plots; plot!(sol2)

alg3 = MethodOfSteps(Tsit5(),constrained=true,max_picard_iters=10,picardabstol=1e-8,picardreltol=1e-10)
sol3 = solve(prob,alg3)
using Plots; plot!(sol3)

alg4 = MethodOfSteps(DP5(),constrained=false,max_picard_iters=100,picardabstol=1e-12,picardreltol=1e-12)
sol4 = solve(prob,alg4)
using Plots; plot(sol4)

@test abs(sol1[end] - sol2[end]) < 1e-3
@test abs(sol1[end] - sol3[end]) < 1e-3
@test abs(sol1[end] - sol4[end]) < 1e-3
