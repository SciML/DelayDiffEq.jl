using DelayDiffEq
using Test

# Test allocation regression for key operations
# These tests verify that in-place variants of critical operations
# remain allocation-free (or have minimal allocations)

@testset "Allocation Regression Tests" begin
    # Setup a simple DDE problem
    function f!(du, u, h, p, t)
        du[1] = -h(p, t - 1.0)[1]
    end

    function h_inplace!(val, p, t; idxs = nothing)
        if idxs === nothing
            val[1] = 1.0
        else
            val[idxs] = 1.0
        end
        return val
    end

    function h(p, t; idxs = nothing)
        if idxs === nothing
            return [1.0]
        else
            return 1.0
        end
    end

    u0 = [1.0]
    tspan = (0.0, 5.0)
    p = nothing

    prob = DDEProblem(f!, u0, h, tspan, p; constant_lags = [1.0])
    alg = MethodOfSteps(Tsit5())

    # Solve to have a dense history
    sol = solve(prob, alg)

    @testset "Solution interpolation (in-place)" begin
        # In-place interpolation should be allocation-free after warmup
        u_cache = zeros(1)

        # Warm up
        sol(u_cache, 2.5)

        # Test for zero allocations (or very minimal)
        allocs = @allocated sol(u_cache, 2.5)
        @test allocs == 0
    end

    @testset "History function (in-place)" begin
        # Create integrator to access history function
        integrator = init(prob, alg)

        hf = integrator.history
        cache_val = zeros(1)

        # Warm up
        hf(cache_val, p, 0.5)

        # In-place history function should be allocation-free
        allocs = @allocated hf(cache_val, p, 0.5)
        @test allocs == 0
    end

    @testset "Solve allocation bounds" begin
        # Test that total solve allocations stay within reasonable bounds
        # This helps detect allocation regressions in the solve path

        # Warm up
        solve(prob, alg)

        # Count allocations for a simple problem
        allocs = @allocated solve(prob, alg)

        # The solve should allocate less than 100KB for this simple problem
        # This is a regression test - if allocations increase significantly,
        # the test will fail
        @test allocs < 100_000

        # Also check number of allocations stays bounded
        # (this is more stable than bytes which depend on array sizes)
        num_allocs = 0
        solve(prob, alg)  # warmup
        stats = @timed solve(prob, alg)
        # We don't have direct access to allocation count in @timed,
        # but we can verify the bytes are bounded
        @test stats.bytes < 100_000
    end

    @testset "Step execution allocation stability" begin
        # For a simple DDE, each step should have bounded allocations
        integrator = init(prob, alg)

        # Take a few steps to warm up
        for _ in 1:5
            step!(integrator)
        end

        # Reset and measure
        reinit!(integrator)
        step!(integrator)  # First step after reinit

        # Measure subsequent step allocations
        allocs_per_step = @allocated step!(integrator)

        # Each step should allocate less than 5KB
        # (most allocations are for growing solution arrays)
        @test allocs_per_step < 5_000
    end
end
