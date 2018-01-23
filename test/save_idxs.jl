using DelayDiffEq, Base.Test

function f_inplace(du,u,h,p,t)
    du[1] = - h(t-1/5)[1] + u[1]
    du[2] = - h(t-1/3)[2] - h(t-1/5)[2]
end
prob_inplace = DDEProblem(f_inplace, t->zeros(2), ones(2), (0.0, 100.0), nothing, [1/5, 1/3])

function f_notinplace(u,h,p,t)
    [-h(t-1/5)[1] + u[1]; -h(t-1/3)[2] - h(t-1/5)[2]]
end
prob_notinplace = DDEProblem(f_notinplace, t->zeros(2), ones(2), (0.0, 100.0), nothing, [1/5, 1/3])

alg = MethodOfSteps(BS3())

for (prob, dense, save_start, save_everystep, saveat) in
    Iterators.product((prob_inplace, prob_notinplace),
                      (true, false),
                      (true, false),
                      (true, false),
                      (Float64[], [25.0, 50.0, 100.0]))

    # save all components (without keyword argument)
    dde_int = init(prob, alg; save_start=save_start, saveat=saveat, dense=dense)
    sol = solve!(dde_int)

    ## solution and solution of ODE integrator contain all components
    @test length(sol.u[1]) == 2
    @test length(dde_int.sol.u[1]) == 2

    ## interpolation
    @test sol(25:100, idxs=2) == [u[1] for u in sol(25:100, idxs=[2])]

    # save all components (with keyword argument)
    dde_int2 = init(prob, alg; save_idxs=[1, 2], save_start=save_start, saveat=saveat,
                    dense=dense)
    sol2 = solve!(dde_int2)

    ## solution and solution of ODE integrator contain all components
    @test length(sol2.u[1]) == 2
    @test length(dde_int2.sol.u[1]) == 2

    ## solution equals solution without keyword arguments
    @test sol.t == sol2.t && sol.u == sol2.u

    ## interpolation
    @test sol(25:100, idxs=2) == sol2(25:100, idxs=2)
    @test sol(25:100, idxs=[2]) == sol2(25:100, idxs=[2])

    # save only second component
    dde_int3 = init(prob, alg; save_idxs=[2], save_start=save_start, saveat=saveat,
                    dense=dense)
    sol3 = solve!(dde_int3)

    ## solution contains only second component
    @test length(sol3.u[1]) == 1

    ## solution of ODE integrator contains both components
    @test length(dde_int3.sol.u[1]) == 2

    ## solution equals second component of complete solution
    @test sol.t == sol3.t && sol[2, :] == sol3[1, :]

    ## interpolation of solution equals second component of interpolation of complete solution
    @test sol(25:100, idxs=2) == sol3(25:100, idxs=1)
    @test sol(25:100, idxs=[2]) == sol3(25:100, idxs=[1])

    # save only second component, scalar index
    dde_int4 = init(prob, alg; save_idxs=2, save_start=save_start, saveat=saveat,
                    dense=dense)
    sol4 = solve!(dde_int4)

    ## solution is only vector of floats
    @test typeof(sol4.u) == Vector{Float64}

    ## solution of ODE integrator contains both components
    @test length(dde_int4.sol.u[1]) == 2

    ## solution equals second component of complete solution
    @test sol.t == sol4.t && sol[2, :] == sol4.u

    ## interpolation of solution equals second component of interpolation of complete solution
    @test sol(25:100, idxs=2) == sol4(25:100, idxs=1)
end
