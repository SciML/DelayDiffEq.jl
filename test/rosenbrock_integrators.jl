using DelayDiffEq, DiffEqProblemLibrary, Base.Test

prob = prob_dde_1delay(1.0)

# ODE algorithms
algs = [Rosenbrock23(), Rosenbrock32(), ROS3P(), Rodas3(),
        RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
        Ros4LStab(), Rodas4(), Rodas42(), Rodas4P(), Rodas5()]

names = ["Rosenbrock23", "Rosenbrock32", "ROS3P", "Rodas3",
         "RosShamp4", "Veldd4", "Velds4", "GRK4T", "GRK4A",
         "Ros4LStab", "Rodas4", "Rodas42", "Rodas4P", "Rodas5"]

for (alg, name) in zip(algs, names)
    print("testing ", name, "... ")
    step_alg = MethodOfSteps(alg)
    solve(prob, step_alg)
    @time solve(prob, step_alg)
end
