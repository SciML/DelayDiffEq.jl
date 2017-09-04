using Unitful, DelayDiffEq, DiffEqProblemLibrary, Base.Test

# Scalar problem, not in-place
prob = DiffEqProblemLibrary.build_prob_dde_1delay_long_scalar_notinplace(1.0u"N", 1.0u"s")

# Unconstrained algorithm without units
alg = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                    fixedpoint_abstol=1e-10, fixedpoint_reltol=1e-4)
sol = solve(prob, alg)

# Unconstrained algorithm with correct units
alg = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                    fixedpoint_abstol=1e-7u"mN", fixedpoint_reltol=1e-4)
sol2 = solve(prob, alg)

@test sol.t == sol2.t && sol.u == sol2.u

# Unconstrained algorithm with incorrect units for absolute tolerance
alg = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                    fixedpoint_abstol=1e-10u"s", fixedpoint_reltol=1e-4)
@test_throws Unitful.DimensionError solve(prob, alg)

# Unconstrained algorithm with units for relative tolerance
alg = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                    fixedpoint_abstol=1e-10, fixedpoint_reltol=1e-4u"N")
@test_throws Unitful.DimensionError solve(prob, alg)

# Constrained algorithm without units
alg = MethodOfSteps(Tsit5(), constrained=true, max_fixedpoint_iters=100,
                    fixedpoint_abstol=1e-10, fixedpoint_reltol=1e-4)
sol = solve(prob, alg)

# Constrained algorithm with correct units
alg = MethodOfSteps(Tsit5(), constrained=true, max_fixedpoint_iters=100,
                    fixedpoint_abstol=1e-7u"mN", fixedpoint_reltol=1e-4)
sol2 = solve(prob, alg)

@test sol.t == sol2.t && sol.u == sol2.u

# Vector problem, in-place
prob = DiffEqProblemLibrary.build_prob_dde_1delay_long(1.0u"N", 1.0u"s")

# Unconstrained algorithm without units
alg = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                    fixedpoint_abstol=1e-10, fixedpoint_reltol=1e-4)
sol = solve(prob, alg)

# Unconstrained algorithm with correct units and both absolute and relative tolerance as
# vector
alg = MethodOfSteps(Tsit5(), constrained=false, max_fixedpoint_iters=100,
                    fixedpoint_abstol=[1e-7u"mN"], fixedpoint_reltol=[1e-4])
sol2 = solve(prob, alg)

@test sol.t == sol2.t && sol.u == sol2.u

# Constrained algorithm without units
alg = MethodOfSteps(Tsit5(), constrained=true, max_fixedpoint_iters=100,
                    fixedpoint_abstol=1e-10, fixedpoint_reltol=1e-4)
sol = solve(prob, alg)

# Constrained algorithm with correct units and both absolute and relative tolerance as
# vector
alg = MethodOfSteps(Tsit5(), constrained=true, max_fixedpoint_iters=100,
                     fixedpoint_abstol=[1e-7u"mN"], fixedpoint_reltol=[1e-4])
sol2 = solve(prob, alg)

@test sol.t == sol2.t && sol.u == sol2.u
