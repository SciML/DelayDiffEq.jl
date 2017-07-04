using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, DataStructures
using Base.Test

@time @testset "Discontinuity Tree Test" begin include("discont_tree_test.jl") end
@time @testset "Constrained Timestep" begin include("constrained.jl") end
@time @testset "Unconstrained Timestep" begin include("unconstrained.jl") end
@time @testset "Events" begin include("events.jl") end
@time @testset "Units" begin include("units.jl") end
@time @testset "Unique Times" begin include("unique_times.jl") end
