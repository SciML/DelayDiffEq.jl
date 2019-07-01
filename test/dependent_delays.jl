include("common.jl")

const alg = MethodOfSteps(BS3())
const dde_int = init(prob_dde_constant_1delay_ip, alg)
const sol = solve!(dde_int)

@testset "reference" begin
  @test sol.errors[:l∞] < 3.0e-5
  @test sol.errors[:final] < 2.1e-5
  @test sol.errors[:l2] < 1.2e-5
end

@testset "constant lag function" begin
  # constant delay specified as lag function
  prob2 = remake(prob_dde_constant_1delay_ip; constant_lags = nothing,
                 dependent_lags = [(u, p, t) -> 1])
  dde_int2 = init(prob2, alg)
  sol2 = solve!(dde_int2)

  @test dde_int.tracked_discontinuities == dde_int2.tracked_discontinuities

  # worse than results above with constant delays specified as scalars
  @test sol2.errors[:l∞] < 4.2e-5
  @test sol2.errors[:final] < 2.2e-5
  @test sol2.errors[:l2] < 1.7e-5

  # simple convergence tests
  @testset "convergence" begin
    sol3 = solve(prob2, alg, abstol=1e-9, reltol=1e-6)

    @test sol3.errors[:l∞] < 7.5e-8
    @test sol3.errors[:final] < 4.6e-8
    @test sol3.errors[:l2] < 3.9e-8

    sol4 = solve(prob2, alg, abstol=1e-13, reltol=1e-13)

    @test sol4.errors[:l∞] < 6.9e-11
    @test sol4.errors[:final] < 1.1e-11
    @test sol4.errors[:l2] < 6.8e-12 # 6.7e-12, relaxed for Win32
  end
end

# without any delays specified is worse
@testset "without delays" begin
  prob2 = remake(prob_dde_constant_1delay_ip; constant_lags = nothing)
  sol2 = solve(prob2, alg)

  @test sol2.errors[:l∞] > 1e-3
  @test sol2.errors[:final] > 4e-5
  @test sol2.errors[:l2] > 7e-4
end
