using DelayDiffEq, OrdinaryDiffEqCore
using DelayDiffEq: DDEVerbosity
using SciMLLogging: SciMLLogging, AbstractMessageLevel
using Test

@testset "DDEVerbosity Construction" begin
    # Test default constructor
    @testset "Default constructor" begin
        v = DDEVerbosity()
        @test v.ode_verbosity isa OrdinaryDiffEqCore.ODEVerbosity
        @test v.discontinuity_tracking isa SciMLLogging.Silent
        @test v.delay_evaluation isa SciMLLogging.Silent
        @test v.constrained_step isa SciMLLogging.Silent
        @test v.residual_control isa SciMLLogging.Silent
        @test v.neutral_delay isa SciMLLogging.Silent
        @test v.state_dependent_delay isa SciMLLogging.Silent
    end

    # Test preset constructors
    @testset "Preset: None" begin
        v = DDEVerbosity(None())
        @test v.ode_verbosity isa OrdinaryDiffEqCore.ODEVerbosity
        @test all(getfield(v, f) isa SciMLLogging.Silent
        for f in fieldnames(DDEVerbosity) if f != :ode_verbosity)
    end

    @testset "Preset: Minimal" begin
        v = DDEVerbosity(Minimal())
        @test v.discontinuity_tracking isa SciMLLogging.Silent
        @test v.delay_evaluation isa SciMLLogging.Silent
        @test v.constrained_step isa SciMLLogging.Silent
        @test v.residual_control isa SciMLLogging.Silent
        @test v.neutral_delay isa SciMLLogging.Silent
        @test v.state_dependent_delay isa SciMLLogging.WarnLevel
    end

    @testset "Preset: Standard" begin
        v = DDEVerbosity(Standard())
        @test v isa DDEVerbosity
    end

    @testset "Preset: Detailed" begin
        v = DDEVerbosity(Detailed())
        @test v.discontinuity_tracking isa SciMLLogging.InfoLevel
        @test v.delay_evaluation isa SciMLLogging.InfoLevel
        @test v.constrained_step isa SciMLLogging.InfoLevel
        @test v.residual_control isa SciMLLogging.InfoLevel
        @test v.neutral_delay isa SciMLLogging.InfoLevel
        @test v.state_dependent_delay isa SciMLLogging.WarnLevel
    end

    @testset "Preset: All" begin
        v = DDEVerbosity(All())
        @test all(getfield(v, f) isa SciMLLogging.InfoLevel
        for f in fieldnames(DDEVerbosity) if f != :ode_verbosity)
    end

    # Test group-level constructor
    @testset "Group-level construction" begin
        v = DDEVerbosity(delay_specific = InfoLevel())
        @test v.discontinuity_tracking isa SciMLLogging.InfoLevel
        @test v.delay_evaluation isa SciMLLogging.InfoLevel
        @test v.constrained_step isa SciMLLogging.InfoLevel
        @test v.residual_control isa SciMLLogging.InfoLevel
        @test v.neutral_delay isa SciMLLogging.InfoLevel
        @test v.state_dependent_delay isa SciMLLogging.InfoLevel
    end

    # Test individual field construction
    @testset "Individual field construction" begin
        v = DDEVerbosity(
            discontinuity_tracking = InfoLevel(),
            delay_evaluation = WarnLevel(),
            state_dependent_delay = DebugLevel()
        )
        @test v.discontinuity_tracking isa SciMLLogging.InfoLevel
        @test v.delay_evaluation isa SciMLLogging.WarnLevel
        @test v.state_dependent_delay isa SciMLLogging.DebugLevel
        @test v.constrained_step isa SciMLLogging.Silent
    end

    # Test mixed group and individual settings
    @testset "Mixed group and individual" begin
        v = DDEVerbosity(
            delay_specific = InfoLevel(),
            state_dependent_delay = WarnLevel()  # Override
        )
        @test v.discontinuity_tracking isa SciMLLogging.InfoLevel
        @test v.delay_evaluation isa SciMLLogging.InfoLevel
        @test v.state_dependent_delay isa SciMLLogging.WarnLevel  # Overridden
    end

    # Test ODE verbosity passthrough
    @testset "ODE verbosity" begin
        ode_v = OrdinaryDiffEqCore.ODEVerbosity(Detailed())
        v = DDEVerbosity(ode_verbosity = ode_v)
        @test v.ode_verbosity === ode_v
    end

    # Test invalid arguments
    @testset "Invalid arguments" begin
        @test_throws ArgumentError DDEVerbosity(invalid_field = InfoLevel())
        @test_throws ArgumentError DDEVerbosity(delay_specific = "not a level")
        @test_throws ArgumentError DDEVerbosity(discontinuity_tracking = 42)
    end
end

@testset "DDEVerbosity Property Access" begin
    # Test direct field access (DDE-specific)
    @testset "Direct DDE field access" begin
        v = DDEVerbosity(discontinuity_tracking = InfoLevel())
        @test v.discontinuity_tracking isa SciMLLogging.InfoLevel
    end

    # Test nested ODE field access via getproperty
    @testset "Nested ODE field access" begin
        ode_v = OrdinaryDiffEqCore.ODEVerbosity(dt_NaN = WarnLevel())
        v = DDEVerbosity(ode_verbosity = ode_v)
        @test v.dt_NaN isa SciMLLogging.WarnLevel
        @test v.max_iters isa SciMLLogging.AbstractMessageLevel
    end

    # Test that invalid field access errors
    @testset "Invalid field access" begin
        v = DDEVerbosity()
        @test_throws ErrorException v.nonexistent_field
    end
end

@testset "DDEVerbosity in solve" begin
    # Define a simple DDE problem
    function dde_func(du, u, h, p, t)
        du[1] = -h(p, t - 1; idxs = 1)
    end
    h(p, t; idxs = nothing) = ones(1)
    prob = DDEProblem(dde_func, [1.0], h, (0.0, 10.0))

    # Test Bool verbose conversion
    @testset "Bool verbose argument" begin
        sol = solve(prob, MethodOfSteps(Tsit5()); verbose = true, dt = 0.1)
        @test sol.retcode == ReturnCode.Success

        sol = solve(prob, MethodOfSteps(Tsit5()); verbose = false, dt = 0.1)
        @test sol.retcode == ReturnCode.Success
    end

    # Test AbstractVerbosityPreset conversion
    @testset "Preset verbose argument" begin
        sol = solve(prob, MethodOfSteps(Tsit5()); verbose = Standard(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success

        sol = solve(prob, MethodOfSteps(Tsit5()); verbose = None(), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
    end

    # Test DDEVerbosity directly
    @testset "DDEVerbosity argument" begin
        v = DDEVerbosity(discontinuity_tracking = InfoLevel())
        @test_logs (:info, r"Handling discontinuity") match_mode=:any solve(
            prob, MethodOfSteps(Tsit5()); verbose = v, dt = 0.1)
    end

    # Test that integrator has correct verbosity type
    @testset "Integrator verbosity type" begin
        integrator = init(prob, MethodOfSteps(Tsit5()); verbose = true, dt = 0.1)
        @test integrator.opts.verbose isa DDEVerbosity

        integrator = init(prob, MethodOfSteps(Tsit5()); verbose = Standard(), dt = 0.1)
        @test integrator.opts.verbose isa DDEVerbosity

        v = DDEVerbosity(delay_evaluation = WarnLevel())
        integrator = init(prob, MethodOfSteps(Tsit5()); verbose = v, dt = 0.1)
        @test integrator.opts.verbose === v
    end
end
