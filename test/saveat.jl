include("common.jl")

const prob = prob_dde_constant_1delay_long_ip
const alg = MethodOfSteps(Tsit5())

# reference integrator and solution
const dde_int = init(prob, alg)
const sol = solve!(dde_int)

@testset "reference" begin
  # solution equals solution of ODE integrator
  @test sol.t == dde_int.sol.t
  @test sol.u == dde_int.sol.u
end

# do not save every step
@testset "not every step (save_start=$save_start)" for save_start in (false, true)
  # for time(s) as scalar (implicitly adds end point as well!) and vectors
  for saveat in (25.0, [25.0, 50.0, 75.0])
    ## minimal ODE solution
    dde_int_min = init(prob, alg; saveat=saveat, save_start=save_start,
                       minimal_solution=true)

    # solution of ODE integrator will be reduced
    @test dde_int_min.saveat !== nothing

    sol_min = solve!(dde_int_min)

    # time point of solution
    @test sol_min.t == (save_start ? [0.0, 25.0, 50.0, 75.0, 100.0] :
                        [25.0, 50.0, 75.0, 100.0])

    # solution of ODE integrator is reduced:
    # [0.0, ≈24.73, ≈25.88, ≈49.43, ≈50.67, ≈74.11, ≈75.34, 100.0]
    @test dde_int_min.sol.t ≈ [0.0, 24.73, 25.88, 49.43, 50.67, 74.11, 75.34, 100.0] atol=1e-2

    # solution lies on interpolation of full solution above
    @test sol(sol_min.t).u == sol_min.u

    ## full ODE solution
    dde_int_full = init(prob, alg; saveat=saveat, save_start=save_start,
                        minimal_solution=false)

    # solution of ODE integrator will not be reduced
    @test dde_int_full.saveat === nothing

    sol_full = solve!(dde_int_full)

    # solution of ODE integrator equals full solution above
    @test sol.t == dde_int_full.sol.t && sol.u == dde_int_full.sol.u

    # solution equals reduced solution above
    @test sol_min.t == sol_full.t && sol_min.u == sol_full.u

    ## dense interpolation
    dde_int_dense = init(prob, alg; saveat=saveat, save_start=save_start,
                         dense=true)

    # solution of ODE integrator will not be reduced
    @test dde_int_dense.saveat === nothing

    sol_dense = solve!(dde_int_dense)

    # time steps of solution
    @test sol_dense.t == (save_start ? [0.0, 25.0, 50.0, 75.0, 100.0] :
                          [25.0, 50.0, 75.0, 100.0])

    # solution of ODE integrator equals full solution above
    @test sol.t == dde_int_dense.sol.t && sol.u == dde_int_dense.sol.u

    # solution lies on interpolation of full solution above
    @test sol(sol_dense.t).u == sol_dense.u

    # full solution above lies on interpolation of solution
    @test sol_dense(sol.t).u == sol.u
  end
end

# save every step
@testset "every step (save_start=$save_start)" for save_start in (false, true)
  # for time(s) as scalar (implicitly adds end point as well!) and vectors
  for saveat in (25.0, [25.0, 50.0, 75.0])
    ## minimal ODE solution (is not possible to minimize ODE solution actually!)
    dde_int_min = init(prob, alg; saveat=saveat, save_everystep=true,
                       save_start=save_start, minimal_solution=true)

    # solution of ODE integrator will not be reduced
    @test dde_int_min.saveat === nothing

    ## full ODE solution
    dde_int_full = init(prob, alg; saveat=saveat, save_everystep=true,
                        save_start=save_start, minimal_solution=false)

    # solution of ODE integrator will not be reduced
    @test dde_int_full.saveat === nothing

    sol_full = solve!(dde_int_full)

    # time steps of solution
    @test symdiff(sol.t, sol_full.t) == (save_start ? [25.0, 50.0, 75.0] :
                                         [0.0, 25.0, 50.0, 75.0])

    # solution of ODE integrator equals full solution above
    @test sol.t == dde_int_full.sol.t && sol.u == dde_int_full.sol.u

    # solution lies on interpolation of full solution above
    @test sol(sol_full.t).u == sol_full.u

    ## dense interpolation
    dde_int_dense = init(prob, alg; saveat=saveat, save_everystep=true,
                         save_start=save_start, dense=true)

    # solution of ODE integrator will not be reduced
    @test dde_int_dense.saveat === nothing

    sol_dense = solve!(dde_int_dense)

    # time steps of solution
    @test symdiff(sol.t, sol_dense.t) == (save_start ? [25.0, 50.0, 75.0] :
                                          [0.0, 25.0, 50.0, 75.0])

    # solution of ODE integrator equals full solution above
    @test sol.t == dde_int_dense.sol.t && sol.u == dde_int_dense.sol.u

    # solution lies on interpolation of full solution above
    @test sol(sol_dense.t).u == sol_dense.u

    # full solution above lies on interpolation of solution
    @test sol_dense(sol.t).u == sol.u
  end
end
