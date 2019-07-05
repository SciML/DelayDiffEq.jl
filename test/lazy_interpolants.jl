include("common.jl")

# simple problems
@testset "simple problems" begin
  prob_ip = prob_dde_constant_1delay_ip
  prob_scalar = prob_dde_constant_1delay_scalar
  ts = 0:0.1:10

  # Vern6
  println("Vern6")
  @testset "Vern6" begin
    alg = MethodOfSteps(Vern6())
    sol_ip = solve(prob_ip, alg)

    @test sol_ip.errors[:l∞] < 8.7e-4
    @test sol_ip.errors[:final] < 7.5e-6
    @test sol_ip.errors[:l2] < 5.5e-4

    sol_scalar = solve(prob_scalar, alg)

    # fails due to floating point issues
    if Sys.WORD_SIZE == 32
      @test_broken sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
      @test_broken sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
    else
      @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)
      @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
    end
  end

  # Vern7
  println("Vern7")
  @testset "Vern7" begin
    alg = MethodOfSteps(Vern7())
    sol_ip = solve(prob_ip, alg)

    @test sol_ip.errors[:l∞] < 4.0e-4
    @test sol_ip.errors[:final] < 3.5e-7
    @test sol_ip.errors[:l2] < 1.9e-4

    sol_scalar = solve(prob_scalar, alg)

    @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)

    # fails due to floating point issues
    if Sys.WORD_SIZE == 32
      @test_broken sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
    else
      @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
    end
  end

  # Vern8
  println("Vern8")
  @testset "Vern8" begin
    alg = MethodOfSteps(Vern8())
    sol_ip = solve(prob_ip, alg)

    @test sol_ip.errors[:l∞] < 2.0e-3
    @test sol_ip.errors[:final] < 1.8e-5
    @test sol_ip.errors[:l2] < 9.0e-4

    sol_scalar = solve(prob_scalar, alg)

    @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)

    # fails due to floating point issues on Win32
    if Sys.WORD_SIZE == 32
      @test_broken sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
    else
      @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
    end
  end

  # Vern9
  println("Vern9")
  @testset "Vern9" begin
    alg = MethodOfSteps(Vern9())
    sol_ip = solve(prob_ip, alg)

    @test sol_ip.errors[:l∞] < 1.5e-3
    @test sol_ip.errors[:final] < 3.8e-6
    @test sol_ip.errors[:l2] < 6.5e-4

    sol_scalar = solve(prob_scalar, alg)

    @test sol_ip(ts, idxs=1) ≈ sol_scalar(ts)

    # fails due to floating point issues
    if Sys.WORD_SIZE == 32
      @test_broken sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
    else
      @test sol_ip.t ≈ sol_scalar.t && sol_ip[1, :] ≈ sol_scalar.u
    end
  end
end

# model of Mackey and Glass
println("Mackey and Glass")
@testset "Mackey and Glass" begin
  prob = prob_dde_DDETST_A1

  # Vern6
  solve(prob, MethodOfSteps(Vern6()))

  # Vern7
  solve(prob, MethodOfSteps(Vern7()))

  # Vern8
  solve(prob, MethodOfSteps(Vern8()))

  # Vern9
  solve(prob, MethodOfSteps(Vern9()))
end
