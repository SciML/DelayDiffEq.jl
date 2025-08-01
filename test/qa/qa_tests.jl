using DelayDiffEq, Test

@testset "Aqua Tests" begin
    using Aqua

    Aqua.test_all(DelayDiffEq; ambiguities = false, piracies = false,
        stale_deps = false, deps_compat = false)
    Aqua.test_ambiguities(DelayDiffEq; recursive = false)
    Aqua.test_stale_deps(DelayDiffEq)
    Aqua.test_deps_compat(DelayDiffEq)
    Aqua.test_piracies(DelayDiffEq)
end

@testset "Explicit Imports Tests" begin
    using ExplicitImports

    @test check_no_implicit_imports(DelayDiffEq; skip = (Base, Core)) === nothing
    @test check_no_stale_explicit_imports(DelayDiffEq) === nothing
    @test check_all_qualified_accesses_via_owners(DelayDiffEq) === nothing
end
