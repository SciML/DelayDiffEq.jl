using DelayDiffEq, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(DelayDiffEq)
    Aqua.test_ambiguities(DelayDiffEq, recursive = false)
    Aqua.test_deps_compat(DelayDiffEq)
    Aqua.test_piracies(DelayDiffEq,
        treat_as_own = [DelayDiffEq.SciMLBase.DESolution,
        DelayDiffEq.SciMLBase.SciMLSolution])
    Aqua.test_project_extras(DelayDiffEq)
    Aqua.test_stale_deps(DelayDiffEq)
    Aqua.test_unbound_args(DelayDiffEq)
    Aqua.test_undefined_exports(DelayDiffEq)
end
