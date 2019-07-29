using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Interface"
  @time @safetestset "Discontinuity Tests" begin include("discontinuities.jl") end
  @time @safetestset "History Function Tests" begin include("history_function.jl") end
  @time @safetestset "Parameterized Function Tests" begin include("parameters.jl") end
  @time @safetestset "Jacobian Tests" begin include("jacobian.jl") end
  @time @safetestset "Return Code Tests" begin include("retcode.jl") end
  @time @safetestset "Composite Solution Tests" begin include("composite_solution.jl") end
  @time @safetestset "Dependent Delay Tests" begin include("dependent_delays.jl") end
  @time @safetestset "Reinitialization Tests" begin include("reinit.jl") end
  @time @safetestset "saveat Tests" begin include("saveat.jl") end
  @time @safetestset "save_idxs Tests" begin include("save_idxs.jl") end
  @time @safetestset "Event Tests" begin include("events.jl") end
  @time @safetestset "Cache Tests" begin include("cache.jl") end
  @time @safetestset "Iterator Tests" begin include("iterator.jl") end
  @time @safetestset "Units Tests" begin include("units.jl") end
  @time @safetestset "Unique Timepoints Tests" begin include("unique_times.jl") end
end

if GROUP == "All" || GROUP == "Integrators"
  @time @safetestset "Constrained Time Steps Tests" begin include("constrained.jl") end
  @time @safetestset "Unconstrained Time Steps Tests" begin include("unconstrained.jl") end
  @time @safetestset "Residual Control Tests" begin include("residual_control.jl") end
  @time @safetestset "Lazy Interpolants Tests" begin include("lazy_interpolants.jl") end
  @time @safetestset "SDIRK Tests" begin include("sdirk_integrators.jl") end
  @time @safetestset "Rosenbrock Tests" begin include("rosenbrock_integrators.jl") end
end
