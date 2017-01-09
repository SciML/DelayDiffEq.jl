using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, DataStructures
using Base.Test

@time @testset "Constrained Timestep" begin include("constrained.jl") end
@time @testset "Unconstrained Timestep" begin include("unconstrained.jl") end
