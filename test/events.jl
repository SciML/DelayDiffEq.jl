using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, Base.Test, DiffEqDevTools, DiffEqCallbacks

lags = [1]
f = function (t,u,h)
    du = - h(t-1)
end
h = (t) -> 0.0

prob = ConstantLagDDEProblem(f, h, 1.0, lags, (0.0, 10.0); iip=DiffEqBase.isinplace(f, 4))
alg = MethodOfSteps(Tsit5(); constrained=false)

# continuous callback

condition = function (t,u,integrator) # Event when event_f(t,u,k) == 0
    t - 2.60
end

affect! = function (integrator)
    integrator.u = -integrator.u
end

cb = ContinuousCallback(condition, affect!)

sol1 = solve(prob, alg, callback=cb)

sol2 = solve(prob, alg, callback=cb, dtmax=0.01)

sol3 = appxtrue(sol1, sol2)

@test sol3.errors[:L2] < 4e-3
@test sol3.errors[:L∞] < 8e-3

# discrete callback

cb = AutoAbstol()

sol1 = solve(prob, alg, callback=cb)

sol2 = solve(prob, alg, callback=cb, dtmax=0.01)

sol3 = appxtrue(sol1, sol2)

@test sol3.errors[:L2] < 3e-2
@test sol3.errors[:L∞] < 7e-2
