using DelayDiffEq, DiffEqProblemLibrary, Base.Test

alg = MethodOfSteps(BS3(); constrained=false)
prob = prob_dde_1delay_scalar_notinplace
integrator = init(prob, alg, dt= 0.01)
solve!(integrator)

u = copy(integrator.sol.u)
t = copy(integrator.sol.t)

reinit!(integrator)
integrator.dt = 0.01
solve!(integrator)

@test u == integrator.sol.u
@test t == integrator.sol.t

integrator = init(prob, alg, dt= 0.01, tstops = [0.5], saveat = [0.33])
sol = solve!(integrator)

u = copy(sol.u)
t = copy(sol.t)

reinit!(integrator)
integrator.dt = 0.01
sol = solve!(integrator)

@test u == sol.u
@test t == sol.t
