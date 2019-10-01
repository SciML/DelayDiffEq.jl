using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Interface"
  @time @safetestset "AD Tests" begin include("interface/ad.jl") end
  @time @safetestset "Backwards Tests" begin include("interface/backwards.jl") end
  @time @safetestset "Composite Solution Tests" begin include("interface/composite_solution.jl") end
  @time @safetestset "Constrained Time Steps Tests" begin include("interface/constrained.jl") end
  @time @safetestset "Dependent Delay Tests" begin include("interface/dependent_delays.jl") end
  @time @safetestset "Discontinuity Tests" begin include("interface/discontinuities.jl") end
  @time @safetestset "Export Tests" begin include("interface/export.jl") end
  @time @safetestset "Fixed-point Iteration Tests" begin include("interface/fpsolve.jl") end
  @time @safetestset "History Function Tests" begin include("interface/history_function.jl") end
  @time @safetestset "Jacobian Tests" begin include("interface/jacobian.jl") end
  @time @safetestset "Parameterized Function Tests" begin include("interface/parameters.jl") end
  @time @safetestset "saveat Tests" begin include("interface/saveat.jl") end
  @time @safetestset "save_idxs Tests" begin include("interface/save_idxs.jl") end
  @time @safetestset "Unconstrained Time Steps Tests" begin include("interface/unconstrained.jl") end
  @time @safetestset "Units Tests" begin include("interface/units.jl") end
end

if GROUP == "All" || GROUP == "Integrators"
  @time @safetestset "Cache Tests" begin include("integrators/cache.jl") end
  @time @safetestset "Event Tests" begin include("integrators/events.jl") end
  @time @safetestset "Iterator Tests" begin include("integrators/iterator.jl") end
  @time @safetestset "Reinitialization Tests" begin include("integrators/reinit.jl") end
  @time @safetestset "Residual Control Tests" begin include("integrators/residual_control.jl") end
  @time @safetestset "Return Code Tests" begin include("integrators/retcode.jl") end
  @time @safetestset "Rosenbrock Tests" begin include("integrators/rosenbrock.jl") end
  @time @safetestset "SDIRK Tests" begin include("integrators/sdirk.jl") end
  @time @safetestset "Unique Timepoints Tests" begin include("integrators/unique_times.jl") end
  @time @safetestset "Verner Tests" begin include("integrators/verner.jl") end
end

if GROUP == "All" || GROUP == "Regression"
  @time @safetestset "Inference Tests" begin include("regression/inference.jl") end
  @time @safetestset "Waltman Problem Tests" begin include("regression/waltman.jl") end
end
