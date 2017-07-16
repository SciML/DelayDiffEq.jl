using DelayDiffEq, OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = prob_dde_1delay_long
alg = MethodOfSteps(Tsit5())

# save at every time step
dde_int = init(prob, alg)
sol = solve!(dde_int)

## solution equals solution of ODE integrator
@test sol.t == dde_int.sol.t
@test sol.u == dde_int.sol.u

# save only at a particular time point
dde_int2 = init(prob, alg; saveat=[50.0])
sol2 = solve!(dde_int2)

## time steps of solution
@test sol2.t == [0.0, 50.0, 100.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int2.sol.t && sol.u == dde_int2.sol.u

## solution lies on interpolation of full solution above
@test sol(sol2.t).u == sol2.u

## full solution above lies on interpolation of solution
@test sol2(sol.t).u == sol.u

# save only at a particular time point and exclude initial time point
dde_int3 = init(prob, alg; saveat=[50.0], save_start=false)
sol3 = solve!(dde_int3)

## time steps of solution
@test sol3.t == [50.0, 100.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int3.sol.t && sol.u == dde_int3.sol.u

## solution lies on interpolation of full solution above
@test sol(sol3.t).u == sol3.u

## full solution above lies on interpolation of solution
@test sol3(sol.t).u == sol.u

# save every step and additionally at a particular time point
dde_int4 = init(prob, alg; saveat=[50.0], save_everystep=true)
sol4 = solve!(dde_int4)

## time steps of solution
@test symdiff(sol.t, sol4.t) == [50.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int4.sol.t && sol.u == dde_int4.sol.u

## solution lies on interpolation of full solution above
@test sol(sol4.t).u == sol4.u

## full solution above lies on interpolation of solution
@test sol4(sol.t).u == sol.u

# save every step, additionally at a particular time point, and exclude initial time point
dde_int5 = init(prob, alg; saveat=[50.0], save_everystep=true, save_start=false)
sol5 = solve!(dde_int5)

## time steps of solution
@test symdiff(sol.t, sol5.t) == [0.0, 50.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int5.sol.t && sol.u == dde_int5.sol.u

## solution lies on interpolation of full solution above
@test sol(sol5.t).u == sol5.u

## full solution above lies on interpolation of solution
@test sol5(sol.t).u == sol.u
