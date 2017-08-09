using DelayDiffEq, Base.Test

f = function (t,u,h,du)
    du[1] = - h(t-1/5)[1] + u[1]
    du[2] = - h(t-1/3)[2] - h(t-1/5)[2]
end
prob = ConstantLagDDEProblem(f, t->zeros(2), ones(2), [1/5, 1/3], (0.0, 100.0))
alg = MethodOfSteps(BS3())

# save all components (without keyword argument)
dde_int = init(prob, alg)
sol = solve!(dde_int)

## solution and solution of ODE integrator contain all components
@test length(sol.u[1]) == 2
@test sol.u == dde_int.sol.u

# save all components (with keyword argument)
dde_int2 = init(prob, alg; save_idxs=[1, 2])
sol2 = solve!(dde_int2)

## solution and solution of ODE integrator contain all components
@test length(sol2.u[1]) == 2
@test sol2.u == dde_int2.sol.u

## solution equals solution without keyword arguments
@test sol.t == sol2.t && sol.u == sol2.u

# save only second component
dde_int3 = init(prob, alg; save_idxs=[2])
sol3 = solve!(dde_int3)

## solution contains only second component
@test length(sol3.u[1]) == 1

## solution equals second component of ODE integrator
@test sol3[1, :] == dde_int3.sol[2, :]

## solution equals second component of complete solution
@test sol.t == sol3.t && sol[2, :] == sol3[1, :]

## interpolation of solution equals second component of interpolation of complete solution
@test sol(0:100, idxs=2) == sol3(0:100, idxs=1)
