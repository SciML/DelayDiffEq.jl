using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, Base.Test

lags = [.2]

f = function (t,u,h,du)
    du[1] = -h(t-.2)[1] + u[1]
end
h = (t) -> [0.0]

prob = ConstantLagDDEProblem(f, h, [1.0], lags, (0.0, 100.0))

for constrained in (false, true)
    alg = MethodOfSteps(Tsit5(), constrained=constrained)
    sol = solve(prob,alg)

    @test allunique(sol.t)
end
