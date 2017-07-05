using Unitful, DelayDiffEq, DiffEqBase, OrdinaryDiffEq, Base.Test

lags = [.2u"s"]

# Scalar problem, not in-place
f = function (t,u,h)
  out = (-h(t-.2u"s") + u) / 1.0u"s"
end
h = (t) -> 0.0u"N"

prob = ConstantLagDDEProblem(f,h,1.0u"N",lags,(0.0u"s",100.0u"s"))

# Unconstrained algorithm without explicit absolut or relative tolerance
alg = MethodOfSteps(Tsit5(),constrained=false,max_fixedpoint_iters=100)
solve(prob,alg)

# Unconstrained algorithm without units
alg = MethodOfSteps(Tsit5(),constrained=false,max_fixedpoint_iters=100,fixedpoint_abstol=1e-12,fixedpoint_reltol=1e-12)
sol = solve(prob,alg)

# Unconstrained algorithm with correct units
alg2 = MethodOfSteps(Tsit5(),constrained=false,max_fixedpoint_iters=100,fixedpoint_abstol=1e-9u"mN",fixedpoint_reltol=1e-12)
sol2 = solve(prob,alg2)

@test sol.t == sol2.t && sol.u == sol2.u

# Unconstrained algorithm with incorrect units for absolute tolerance
alg = MethodOfSteps(Tsit5(),constrained=false,max_fixedpoint_iters=100,fixedpoint_abstol=1e-12u"s",fixedpoint_reltol=1e-12)
@test_throws Unitful.DimensionError solve(prob,alg)

# Unconstrained algorithm with units for relative tolerance
alg = MethodOfSteps(Tsit5(),constrained=false,max_fixedpoint_iters=100,fixedpoint_abstol=1e-12,fixedpoint_reltol=1e-12u"N")
@test_throws Unitful.DimensionError solve(prob,alg)

# Constrained algorithm without explicit absolute or relative tolerance
alg = MethodOfSteps(Tsit5(),constrained=true,max_fixedpoint_iters=100)
solve(prob,alg)

# Constrained algorithm without units
alg = MethodOfSteps(Tsit5(),constrained=true,max_fixedpoint_iters=100,fixedpoint_abstol=1e-12,fixedpoint_reltol=1e-12)
sol = solve(prob,alg)

# Constrained algorithm with correct units
alg2 = MethodOfSteps(Tsit5(),constrained=true,max_fixedpoint_iters=100,fixedpoint_abstol=1e-9u"mN",fixedpoint_reltol=1e-12)
sol2 = solve(prob,alg2)

@test sol.t == sol2.t && sol.u == sol2.u

# Vector problem, in-place
f = function (t,u,h,du)
  du[1] = (-h(t-.2u"s")[1] + u[1]) / 1.0u"s"
end
h = (t) -> [0.0u"N"]

prob = ConstantLagDDEProblem(f,h,[1.0u"N"],lags,(0.0u"s",100.0u"s"))

# Unconstrained algorithm without explicit absolute or relative tolerance
alg = MethodOfSteps(Tsit5(),constrained=false,max_fixedpoint_iters=100)
solve(prob,alg)

# Unconstrained algorithm without units
alg = MethodOfSteps(Tsit5(),constrained=false,max_fixedpoint_iters=100,fixedpoint_abstol=1e-12,fixedpoint_reltol=1e-12)
sol = solve(prob,alg)

# Unconstrained algorithm with correct units and both absolute and relative tolerance as vector
alg2 = MethodOfSteps(Tsit5(),constrained=false,max_fixedpoint_iters=100,fixedpoint_abstol=[1e-9u"mN"],fixedpoint_reltol=[1e-12])
sol2 = solve(prob,alg2)

@test sol.t == sol2.t && sol.u == sol2.u

# Constrained algorithm without explicit absolute or relative tolerance
alg = MethodOfSteps(Tsit5(),constrained=true,max_fixedpoint_iters=100)
solve(prob,alg)

# Constrained algorithm without units
alg = MethodOfSteps(Tsit5(),constrained=true,max_fixedpoint_iters=100,fixedpoint_abstol=1e-12,fixedpoint_reltol=1e-12)
sol = solve(prob,alg)

# Constrained algorithm with correct units and both absolute and relative tolerance as vector
alg2 = MethodOfSteps(Tsit5(),constrained=true,max_fixedpoint_iters=100,fixedpoint_abstol=[1e-9u"mN"],fixedpoint_reltol=[1e-12])
sol2 = solve(prob,alg2)

@test sol.t == sol2.t && sol.u == sol2.u
