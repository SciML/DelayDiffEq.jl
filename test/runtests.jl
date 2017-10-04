using Base.Test

@time @testset "Discontinuity Tests" begin include("discontinuities.jl") end
@time @testset "HistoryFunction Tests" begin include("history_function.jl") end
@time @testset "Constrained Timestep" begin include("constrained.jl") end
@time @testset "Unconstrained Timestep" begin include("unconstrained.jl") end
@time @testset "Dependent Delay Tests" begin include("dependent_delays.jl") end
@time @testset "Saveat Test" begin include("saveat.jl") end
@time @testset "Save_idxs Test" begin include("save_idxs.jl") end
@time @testset "Events" begin include("events.jl") end
@time @testset "Units" begin include("units.jl") end
@time @testset "Unique Times" begin include("unique_times.jl") end
@time @testset "Residual Control" begin include("residual_control.jl") end
@time @testset "Lazy Interpolants" begin include("lazy_interpolants.jl") end
@time @testset "SDIRK Integrators" begin include("sdirk_integrators.jl") end
@time @testset "Rosenbrock Integrators" begin include("rosenbrock_integrators.jl") end
