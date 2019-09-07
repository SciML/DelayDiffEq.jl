using DelayDiffEq
using Test

@testset "in-place" begin
  # define functions (Hutchinson's equation)
  function f(du, u, h, p, t)
    du[1] = u[1] * (1  - h(p, t - 1)[1])
    nothing
  end

  njacs = Ref(0)
  function jac(J, u, h, p, t)
    njacs[] += 1
    J[1, 1] = 1 - h(p, t - 1)[1]
    nothing
  end

  nWfacts = Ref(0)
  function Wfact(W, u, h, p, dtgamma, t)
    nWfacts[] += 1
    W[1,1] = dtgamma * (1 - h(p, t - 1)[1]) - 1
    nothing
  end

  nWfact_ts = Ref(0)
  function Wfact_t(W, u, h, p, dtgamma, t)
    nWfact_ts[] += 1
    W[1,1] = 1 - h(p, t - 1)[1] - inv(dtgamma)
    nothing
  end

  h(p, t) = [0.0]

  # define problems
  prob = DDEProblem(DDEFunction{true}(f), [1.0], h, (0.0, 40.0); constant_lags = [1])
  prob_jac = remake(prob; f = DDEFunction{true}(f; jac = jac))
  prob_Wfact = remake(prob; f = DDEFunction{true}(f; Wfact = Wfact))
  prob_Wfact_t = remake(prob; f = DDEFunction{true}(f; Wfact_t = Wfact_t))

  # compute solutions
  for alg in (Rosenbrock23(), TRBDF2())
    sol = solve(prob, MethodOfSteps(alg))

    ## Jacobian
    njacs[] = 0
    sol_jac = solve(prob_jac, MethodOfSteps(alg))
    
    # check number of function evaluations
    @test !iszero(njacs[])
    @test njacs[] == sol_jac.destats.njacs
    if alg isa Rosenbrock23
      @test njacs[] == sol_jac.destats.nw
    else
      @test_broken njacs[] == sol_jac.destats.nw
    end

    # check resulting solution
    @test sol.t ≈ sol_jac.t
    @test sol.u ≈ sol_jac.u

    ## Wfact
    nWfacts[] = 0
    sol_Wfact = solve(prob_Wfact, MethodOfSteps(alg))

    # check number of function evaluations
    if alg isa Rosenbrock23
      @test !iszero(nWfacts[])
      @test nWfacts[] == njacs[]
      @test iszero(sol_Wfact.destats.njacs)
    else
      @test_broken !iszero(nWfacts[])
      @test_broken nWfacts[] == njacs[]
      @test_broken iszero(sol_Wfact.destats.njacs)
    end
    @test_broken nWfacts[] == sol_Wfact.destats.nw

    # check resulting solution
    @test sol.t ≈ sol_Wfact.t
    @test sol.u ≈ sol_Wfact.u

    ## Wfact_t
    nWfact_ts[] = 0
    sol_Wfact_t = solve(prob_Wfact_t, MethodOfSteps(alg))

    # check number of function evaluations
    if alg isa Rosenbrock23
      @test_broken !iszero(nWfact_ts[])
      @test_broken nWfact_ts[] == njacs[]
      @test_broken iszero(sol_Wfact_t.destats.njacs)
    else
      @test !iszero(nWfact_ts[])
      @test_broken nWfact_ts[] == njacs[]
      @test iszero(sol_Wfact_t.destats.njacs)
    end
    @test_broken nWfact_ts[] == sol_Wfact_t.destats.nw

    # check resulting solution
    if alg isa Rosenbrock23
      @test sol.t ≈ sol_Wfact_t.t
      @test sol.u ≈ sol_Wfact_t.u
    else
      @test_broken sol.t ≈ sol_Wfact_t.t
      @test_broken sol.u ≈ sol_Wfact_t.u
    end
  end
end

@testset "out-of-place" begin
  # define functions (Hutchinson's equation)
  f(u, h, p, t) = u[1] .* (1  .- h(p, t - 1))

  njacs = Ref(0)
  function jac(u, h, p, t)
    njacs[] += 1
    reshape(1 .- h(p, t - 1), 1, 1)
  end

  nWfacts = Ref(0)
  function Wfact(u, h, p, dtgamma, t)
    nWfacts[] += 1
    reshape(dtgamma .* (1 .- h(p, t - 1)) .- 1, 1, 1)
  end

  nWfact_ts = Ref(0)
  function Wfact_t(u, h, p, dtgamma, t)
    nWfact_ts[] += 1
    reshape((1 - inv(dtgamma)) .- h(p, t - 1), 1, 1)
  end

  h(p, t) = [0.0]

  # define problems
  prob = DDEProblem(DDEFunction{false}(f), [1.0], h, (0.0, 40.0); constant_lags = [1])
  prob_jac = remake(prob; f = DDEFunction{false}(f; jac = jac))
  prob_Wfact = remake(prob; f = DDEFunction{false}(f; Wfact = Wfact))
  prob_Wfact_t = remake(prob; f = DDEFunction{false}(f; Wfact_t = Wfact_t))

  # compute solutions
  for alg in (Rosenbrock23(), TRBDF2())
    sol = solve(prob, MethodOfSteps(alg))
    
    ## Jacobian
    njacs[] = 0
    sol_jac = solve(prob_jac, MethodOfSteps(alg))

    # check number of function evaluations
    @test !iszero(njacs[])
    @test_broken njacs[] == sol_jac.destats.njacs
    @test_broken njacs[] == sol_jac.destats.nw

    # check resulting solution
    @test sol.t ≈ sol_jac.t
    @test sol.u ≈ sol_jac.u

    ## Wfact
    nWfacts[] = 0
    sol_Wfact = solve(prob_Wfact, MethodOfSteps(alg))

    # check number of function evaluations
    @test_broken !iszero(nWfacts[])
    @test_broken nWfacts[] == njacs[]
    @test_broken iszero(sol_Wfact.destats.njacs)
    @test_broken nWfacts[] == sol_Wfact.destats.nw

    # check resulting solution
    @test sol.t ≈ sol_Wfact.t
    @test sol.u ≈ sol_Wfact.u

    ## Wfact_t
    nWfact_ts[] = 0
    sol_Wfact_t = solve(prob_Wfact_t, MethodOfSteps(alg))

    # check number of function evaluations
    @test_broken !iszero(nWfact_ts[])
    @test_broken nWfact_ts[] == njacs[]
    @test_broken iszero(sol_Wfact_ts.destats.njacs)
    @test_broken nWfact_ts[] == sol_Wfact_t.destats.nw

    # check resulting solution
    @test sol.t ≈ sol_Wfact_t.t
    @test sol.u ≈ sol_Wfact_t.u
  end
end
