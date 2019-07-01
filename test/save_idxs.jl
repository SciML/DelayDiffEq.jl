include("common.jl")

# out-of-place problem
function f_notinplace(u,h,p,t)
  [-h(p, t-1/5)[1] + u[1]; -h(p, t-1/3)[2] - h(p, t-1/5)[2]]
end
const prob_notinplace = DDEProblem(f_notinplace, ones(2), (p, t)->zeros(2), (0.0, 100.0),
                                   constant_lags = [1/5, 1/3])

# in-place problem
function f_inplace(du,u,h,p,t)
  du[1] = - h(p, t-1/5)[1] + u[1]
  du[2] = - h(p, t-1/3)[2] - h(p, t-1/5)[2]
end
const prob_inplace = DDEProblem(f_inplace, ones(2), (p, t)->zeros(2), (0.0, 100.0),
                                constant_lags = [1/5, 1/3])

const alg = MethodOfSteps(BS3())

@testset for prob in (prob_notinplace, prob_inplace),
  dense in (true, false), save_start in (true, false),
  save_everystep in (true, false),
  saveat in (Float64[], [25.0, 50.0, 100.0])

  # reference solution
  dde_int = init(prob, alg; save_start=save_start, saveat=saveat,
                 dense=dense)
  sol = solve!(dde_int)

  # save all components
  @testset "all components" begin
    # without keyword argument
    @testset "without keyword" begin
      ## solution and solution of ODE integrator contain all components
      @test length(sol.u[1]) == 2
      @test length(dde_int.sol.u[1]) == 2

      ## interpolation
      @test sol(25:100, idxs=2) ≈ [u[1] for u in sol(25:100, idxs=[2])]
    end

    # with keyword argument
    @testset "with keyword" begin
      dde_int2 = init(prob, alg; save_idxs=[1, 2], save_start=save_start,
                      saveat=saveat, dense=dense)
      sol2 = solve!(dde_int2)

      ## solution and solution of ODE integrator contain all components
      @test length(sol2.u[1]) == 2
      @test length(dde_int2.sol.u[1]) == 2

      ## solution equals solution without keyword arguments
      @test sol.t == sol2.t && sol.u == sol2.u

      ## interpolation
      @test sol(25:100, idxs=2) ≈ sol2(25:100, idxs=2)
      @test sol(25:100, idxs=[2]) ≈ sol2(25:100, idxs=[2])
    end
  end

  # save only second component
  @testset "second component" begin
    # array index
    @testset "array index" begin
      dde_int2 = init(prob, alg; save_idxs=[2], save_start=save_start,
                      saveat=saveat, dense=dense)
      sol2 = solve!(dde_int2)

      ## solution contains only second component
      @test length(sol2.u[1]) == 1

      ## solution of ODE integrator contains both components
      @test length(dde_int2.sol.u[1]) == 2

      ## solution equals second component of complete solution
      @test sol.t == sol2.t && sol[2, :] == sol2[1, :]

      ## interpolation of solution equals second component of
      ## interpolation of complete solution
      @test sol(25:100, idxs=2) ≈ sol2(25:100, idxs=1)
      @test sol(25:100, idxs=[2]) ≈ sol2(25:100, idxs=[1])
    end

    # scalar index
    @testset "scalar index" begin
      dde_int2 = init(prob, alg; save_idxs=2, save_start=save_start,
                      saveat=saveat, dense=dense)
      sol2 = solve!(dde_int2)

      ## solution is only vector of floats
      @test typeof(sol2.u) == Vector{Float64}

      ## solution of ODE integrator contains both components
      @test length(dde_int2.sol.u[1]) == 2

      ## solution equals second component of complete solution
      @test sol.t ≈ sol2.t && sol[2, :] ≈ sol2.u

      ## interpolation of solution equals second component of
      ## interpolation of complete solution
      @test sol(25:100, idxs=2) ≈ sol2(25:100, idxs=1)
    end
  end
end
