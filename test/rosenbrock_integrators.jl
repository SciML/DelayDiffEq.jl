include("common.jl")
@testset "Rosenbrock integrators" begin
    prob_inplace = prob_dde_1delay
    prob_notinplace = prob_dde_1delay_scalar_notinplace

    # ODE algorithms
    algs = [Rosenbrock23(), Rosenbrock32(), ROS3P(), Rodas3(),
            RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
            Ros4LStab(), Rodas4(), Rodas42(), Rodas4P(), Rodas5()]

    @testset for alg in algs
        @show alg
        stepsalg = MethodOfSteps(alg)
        solve(prob_inplace, stepsalg)
        solve(prob_notinplace, stepsalg)
    end
end
