using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, DataStructures
using Base.Test

@time @testset "Constrained Timestep" begin include("constrained.jl") end
