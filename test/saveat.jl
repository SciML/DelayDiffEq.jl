using DelayDiffEq, OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = prob_dde_1delay_long
alg = MethodOfSteps(Tsit5())

# save at every time step
dde_int = init(prob, alg)
sol = solve!(dde_int)

## full solution equals both solution of DDE integrator and of ODE integrator
@test sol.t == dde_int.sol.t && sol.u == dde_int.sol.u
@test sol.t == dde_int.integrator.sol.t && sol.u == dde_int.integrator.sol.u

# save only at a particular time point
dde_int2 = init(prob, alg; saveat=[50.0])
sol2 = solve!(dde_int2)

## time steps of partial solution
@test sol2.t == [0.0, 50.0, 100.0]

## partial solution equals solution of DDE integrator
@test sol2.t == dde_int2.sol.t && sol2.u == dde_int2.sol.u

## solution of ODE integrator equals full solution
@test sol.t == dde_int2.integrator.sol.t && sol.u == dde_int2.integrator.sol.u

## partial solution equals interpolation of full solution
@test sol(sol2.t).u == sol2.u

# save only at a particular time point and exclude initial time point
dde_int3 = init(prob, alg; saveat=[50.0], save_start=false)
sol3 = solve!(dde_int3)

## time steps of partial solution
@test sol3.t == [50.0, 100.0]

## partial solution equals solution of DDE integrator
@test sol3.t == dde_int3.sol.t && sol3.u == dde_int3.sol.u

## solution of ODE integrator equals full solution
@test sol.t == dde_int3.integrator.sol.t && sol.u == dde_int3.integrator.sol.u

## partial solution equals interpolation of full solution
@test sol(sol3.t).u == sol3.u
