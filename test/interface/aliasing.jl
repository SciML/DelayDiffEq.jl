using DelayDiffEq, DDEProblemLibrary

# For now, testing if the old keyword works the same as the new alias keyword
prob_ip = prob_dde_constant_1delay_ip
prob_scalar = prob_dde_constant_1delay_scalar
ts = 0:0.1:10

noreuse = NLNewton(fast_convergence_cutoff = 0)

const working_algs = [ImplicitMidpoint(), SSPSDIRK2(), KenCarp5(nlsolve = noreuse),
    ImplicitEuler(nlsolve = noreuse), Trapezoid(nlsolve = noreuse),
    TRBDF2(nlsolve = noreuse)]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in working_algs
    println(nameof(typeof(alg)))

    stepsalg = MethodOfSteps(alg)
    sol_new_alias = solve(prob_ip, stepsalg; dt = 0.1, alias = DDEAliasSpecifier(alias_u0 = true))
    sol_old_alias = solve(
        prob_ip, stepsalg; dt = 0.1, alias = alias_u0 = true)

    @test sol_new_alias == sol_old_alias
end