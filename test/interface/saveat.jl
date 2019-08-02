using DelayDiffEq, DiffEqProblemLibrary.DDEProblemLibrary
using Test

DDEProblemLibrary.importddeproblems()

const prob = DDEProblemLibrary.prob_dde_constant_1delay_long_ip
const alg = MethodOfSteps(Tsit5())

# reference integrator and solution
const dde_int = init(prob, alg)
const sol = solve!(dde_int)

@testset "reference" begin
  # solution equals solution of DDE integrator
  @test sol.t == dde_int.sol.t
  @test sol.u == dde_int.sol.u

  # solution equals solution of ODE integrator
  @test sol.t == dde_int.integrator.sol.t
  @test sol.u == dde_int.integrator.sol.u
end

# do not save every step
@testset "not every step (save_start=$save_start)" for save_start in (false, true)
  # for time(s) as scalar (implicitly adds end point as well!) and vectors
  for saveat in (25.0, [25.0, 50.0, 75.0])
    dde_int2 = init(prob, alg; saveat=saveat, save_start=save_start)

    # end point is saved if saveat is a scalar
    @test dde_int2.opts.save_end == (saveat isa Number)

    sol2 = solve!(dde_int2)

    # solution is equal to solution of DDE integrator
    @test sol2.t == dde_int2.sol.t
    @test sol2.u == dde_int2.sol.u

    # time point of solution
    if saveat isa Number
      @test sol2.t == (save_start ? [0.0, 25.0, 50.0, 75.0, 100.0] : [25.0, 50.0, 75.0, 100.0])
    else
      @test sol2.t == (save_start ? [0.0, 25.0, 50.0, 75.0] : [25.0, 50.0, 75.0])
    end

    # history is equal to solution above
    @test sol.t == dde_int2.integrator.sol.t
    @test sol.u == dde_int2.integrator.sol.u
  end
end

# do not save every step
@testset "not every step (save_end=$save_end)" for save_end in (false, true)
  # for time(s) as scalar (implicitly adds end point as well!) and vectors
  for saveat in (25.0, [25.0, 50.0, 75.0])
    dde_int2 = init(prob, alg; saveat=saveat, save_end=save_end)

    # start point is saved if saveat is a scalar
    @test dde_int2.opts.save_start == (saveat isa Number)

    sol2 = solve!(dde_int2)

    # solution is equal to solution of DDE integrator
    @test sol2.t == dde_int2.sol.t
    @test sol2.u == dde_int2.sol.u

    # time point of solution
    if saveat isa Number
      @test sol2.t == [0.0, 25.0, 50.0, 75.0, 100.0]
    else
      @test sol2.t == (save_end ? [25.0, 50.0, 75.0, 100.0] : [25.0, 50.0, 75.0])
    end

    # history is equal to solution above
    @test sol.t == dde_int2.integrator.sol.t
    @test sol.u == dde_int2.integrator.sol.u
  end
end

# save every step
@testset "every step (save_start=$save_start)" for save_start in (false, true)
  for saveat in (25.0, [25.0, 50.0, 75.0])
    dde_int2 = init(prob, alg; saveat=saveat, save_everystep=true,
                    save_start=save_start)

    # end point is saved implicitly
    @test dde_int2.opts.save_end

    sol2 = solve!(dde_int2)

    # solution is equal to solution of DDE integrator
    @test sol2.t == dde_int2.sol.t
    @test sol2.u == dde_int2.sol.u

    # time points of solution
    if saveat isa Number
      @test symdiff(sol.t, sol2.t) == (save_start ? [100.0, 25.0, 50.0, 75.0] :
                                       [0.0, 100.0, 25.0, 50.0, 75.0])
    else
      @test symdiff(sol.t, sol2.t) == (save_start ? [25.0, 50.0, 75.0] :
                                       [0.0, 25.0, 50.0, 75.0])
    end

    # history is equal to solution above
    @test sol.t == dde_int2.integrator.sol.t
    @test sol.u == dde_int2.integrator.sol.u
  end
end

# save every step
@testset "every step (save_end=$save_end)" for save_end in (false, true)
  for saveat in (25.0, [25.0, 50.0, 75.0])
    dde_int2 = init(prob, alg; saveat=saveat, save_everystep=true,
                    save_end=save_end)

    # start point is saved implicitly
    @test dde_int2.opts.save_start

    sol2 = solve!(dde_int2)

    # solution is equal to solution of DDE integrator
    @test sol2.t == dde_int2.sol.t
    @test sol2.u == dde_int2.sol.u

    # time points of solution
    if saveat isa Number
      @test symdiff(sol.t, sol2.t) == (save_end ? [100.0, 25.0, 50.0, 75.0] :
                                       [100.0, 25.0, 50.0, 75.0])
    else
      @test symdiff(sol.t, sol2.t) == (save_end ? [25.0, 50.0, 75.0] :
                                       [25.0, 50.0, 75.0])
    end

    # history is equal to solution above
    @test sol.t == dde_int2.integrator.sol.t
    @test sol.u == dde_int2.integrator.sol.u
  end
end
