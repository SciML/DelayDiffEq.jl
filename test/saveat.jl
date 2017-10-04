using DelayDiffEq, DiffEqProblemLibrary, Base.Test

prob = prob_dde_1delay_long
alg = MethodOfSteps(Tsit5())

# save at every time step
dde_int = init(prob, alg)
sol = solve!(dde_int)

## solution equals solution of ODE integrator
@test sol.t == dde_int.sol.t
@test sol.u == dde_int.sol.u

# save only at a particular time point
dde_int2 = init(prob, alg; saveat=[25.0, 50.0, 75.0])

## solution of ODE integrator will be reduced
@test dde_int2.saveat != nothing

sol2 = solve!(dde_int2)

## time steps of solution
@test sol2.t == [0.0, 25.0, 50.0, 75.0, 100.0]

## solution of ODE integrator is reduced:
## [0.0, ≈24.67, ≈25.89, ≈49.23, ≈50.47, ≈73.91, ≈75.14, 100.0]
@test dde_int2.sol.t ≈ [0.0, 24.67, 25.89, 49.23, 50.47, 73.91, 75.14, 100.0] atol=1.1e-2

## solution lies on interpolation of full solution above
@test sol(sol2.t).u == sol2.u

# save only at a particular time point, force full ODE solution
dde_int2_full = init(prob, alg; saveat=[25.0, 50.0, 75.0], minimal_solution=false)

## solution of ODE integrator will not be reduced
@test dde_int2_full.saveat == nothing

sol2_full = solve!(dde_int2_full)

## solution of ODE integrator equals full solution above
@test sol.t == dde_int2_full.sol.t && sol.u == dde_int2_full.sol.u

## solution equals reduced solution above
@test sol2.t == sol2_full.t && sol2.u == sol2_full.u

# save only at a particular time point (dense interpolation)
dde_int2_dense = init(prob, alg; saveat=[25.0, 50.0, 75.0], dense=true)

## solution of ODE integrator will not be reduced
@test dde_int2_dense.saveat == nothing

sol2_dense = solve!(dde_int2_dense)

## time steps of solution
@test sol2_dense.t == [0.0, 25.0, 50.0, 75.0, 100.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int2_dense.sol.t && sol.u == dde_int2_dense.sol.u

## solution lies on interpolation of full solution above
@test sol(sol2_dense.t).u == sol2_dense.u

## full solution above lies on interpolation of solution
@test sol2_dense(sol.t).u == sol.u

# save only at a particular time point and exclude initial time point
dde_int3 = init(prob, alg; saveat=[25.0, 50.0, 75.0], save_start=false)

## solution of ODE integrator will be reduced
@test dde_int3.saveat != nothing

sol3 = solve!(dde_int3)

## time steps of solution
@test sol3.t == [25.0, 50.0, 75.0, 100.0]

## solution of ODE integrator equals reduced solution of ODE integrator above
@test dde_int3.sol.t == dde_int2.sol.t

## solution lies on interpolation of full solution above
@test sol(sol3.t).u == sol3.u

# save only at a particular time point and exclude initial time point (dense interpolation)
dde_int3_dense = init(prob, alg; saveat=[25.0, 50.0, 75.0], save_start=false, dense=true)
sol3_dense = solve!(dde_int3_dense)

## time steps of solution
@test sol3_dense.t == [25.0, 50.0, 75.0, 100.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int3_dense.sol.t && sol.u == dde_int3_dense.sol.u

## solution lies on interpolation of full solution above
@test sol(sol3_dense.t).u == sol3_dense.u

## full solution above lies on interpolation of solution
@test sol3_dense(sol.t).u == sol.u

# save every step and additionally at a particular time point
dde_int4 = init(prob, alg; saveat=[25.0, 50.0, 75.0], save_everystep=true)
sol4 = solve!(dde_int4)

## time steps of solution
@test symdiff(sol.t, sol4.t) == [25.0, 50.0, 75.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int4.sol.t && sol.u == dde_int4.sol.u

## solution lies on interpolation of full solution above
@test sol(sol4.t).u == sol4.u

# save every step and additionally at a particular time point (dense interpolation)
dde_int4_dense = init(prob, alg; saveat=[25.0, 50.0, 75.0], save_everystep=true, dense=true)
sol4_dense = solve!(dde_int4_dense)

## time steps of solution
@test symdiff(sol.t, sol4_dense.t) == [25.0, 50.0, 75.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int4_dense.sol.t && sol.u == dde_int4_dense.sol.u

## solution lies on interpolation of full solution above
@test sol(sol4_dense.t).u == sol4_dense.u

## full solution above lies on interpolation of solution
@test sol4_dense(sol.t).u == sol.u

# save every step, additionally at a particular time point, and exclude initial time point
dde_int5 = init(prob, alg; saveat=[25.0, 50.0, 75.0], save_everystep=true, save_start=false)
sol5 = solve!(dde_int5)

## time steps of solution
@test symdiff(sol.t, sol5.t) == [0.0, 25.0, 50.0, 75.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int5.sol.t && sol.u == dde_int5.sol.u

## solution lies on interpolation of full solution above
@test sol(sol5.t).u == sol5.u

# save every step, additionally at a particular time point, and exclude initial time point
# (dense interpolation)
dde_int5_dense = init(prob, alg; saveat=[25.0, 50.0, 75.0], save_everystep=true,
                      save_start=false, dense=true)
sol5_dense = solve!(dde_int5_dense)

## time steps of solution
@test symdiff(sol.t, sol5_dense.t) == [0.0, 25.0, 50.0, 75.0]

## solution of ODE integrator equals full solution above
@test sol.t == dde_int5_dense.sol.t && sol.u == dde_int5_dense.sol.u

## solution lies on interpolation of full solution above
@test sol(sol5_dense.t).u == sol5_dense.u

## full solution above lies on interpolation of solution
@test sol5_dense(sol.t).u == sol.u
