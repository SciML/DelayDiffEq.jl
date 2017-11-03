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
