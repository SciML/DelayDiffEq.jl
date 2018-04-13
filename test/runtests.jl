using DelayDiffEq, DiffEqProblemLibrary, Base.Test

tests = ["discontinuities.jl",
         "history_function.jl",
         "constrained.jl",
         "unconstrained.jl",
         "parameters.jl",
         "retcode.jl",
         "composite_solution.jl",
         "dependent_delays.jl",
         "reinit.jl",
         "saveat.jl",
         "save_idxs.jl",
         "events.jl",
         "units.jl",
         "unique_times.jl",
         "residual_control.jl",
         "lazy_interpolants.jl",
         "sdirk_integrators.jl",
         "rosenbrock_integrators.jl"]

for test in tests
    @time include(test)
end
