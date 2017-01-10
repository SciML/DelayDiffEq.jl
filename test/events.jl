using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, Base.Test, DiffEqDevTools

lags = [1]
f = function (t,u,h)
  du = - h(t-1)
end
h = (t) -> 0.0


prob = DDEProblem(f,h,1.0,lags,(0.0,10.0);iip=DiffEqBase.isinplace(f,4))
alg = MethodOfSteps(Tsit5();constrained=false)

condtion= function (t,u,integrator) # Event when event_f(t,u,k) == 0
  t - 2.60
end

affect! = function (integrator)
  integrator.u = -integrator.u
end

rootfind = true
save_positions = (true,true)
callback = Callback(condtion,affect!,rootfind,save_positions)

sol = solve(prob,alg,callback=callback)

sol2= solve(prob,alg,callback=callback,dtmax=0.01)

sol3 = appxtrue(sol,sol2)

@test sol3.errors[:L2] < 4e-3
@test sol3.errors[:L∞] < 8e-3
