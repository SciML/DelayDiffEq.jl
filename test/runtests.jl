using DelayDiffEqPaper

using CSV
using DataFrames
using DDEProblemLibrary
using DelayDiffEq
using DiffEqBase
using DiffEqDevTools
using FiniteDiff

using Statistics
using Test

@testset "DelayDiffEqPaper" begin
    @testset "Anderson" begin include("anderson.jl") end
    @testset "Benchmark" begin include("benchmark.jl") end
    @testset "Mackey-Glass" begin include("mackey_glass.jl") end
    @testset "Sensitivity" begin include("sensitivity.jl") end
    @testset "Waltman" begin include("waltman.jl") end
end
