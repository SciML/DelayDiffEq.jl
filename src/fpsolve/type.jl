# solver

abstract type AbstractFPSolver <: DiffEqBase.AbstractNLSolver end

mutable struct FPSolver{algType<:DiffEqBase.AbstractNLSolverAlgorithm,IIP,uTolType,C} <: AbstractFPSolver
  alg::algType
  κ::uTolType
  η::uTolType
  ndz::uTolType
  fast_convergence_cutoff::uTolType
  iter::Int
  maxiters::Int
  status::DiffEqBase.NLStatus
  cache::C
end

# caches

abstract type AbstractFPSolverCache <: DiffEqBase.AbstractNLSolverCache end

struct FPFunctionalCache{uNoUnitsType} <: AbstractFPSolverCache
  atmp::uNoUnitsType
end

struct FPFunctionalConstantCache <: AbstractFPSolverCache end

mutable struct FPAndersonCache{uType,uNoUnitsType,uEltypeNoUnits,D} <: AbstractFPSolverCache
  """residuals `g(z) - z` of fixed-point iteration"""
  dz::uType
  atmp::uNoUnitsType
  """value `g(z)` of previous fixed-point iteration"""
  gzprev::uType
  """residuals `g(z) - z` of previous fixed-point iteration"""
  dzprev::uType
  Δgzs::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  droptol::D
end

mutable struct FPAndersonConstantCache{uType,uEltypeNoUnits,D} <: AbstractFPSolverCache
  """residuals `g(z) - z` of fixed-point iteration"""
  dz::uType
  """value `g(z)` of previous fixed-point iteration"""
  gzprev::uType
  """residuals `g(z) - z` of previous fixed-point iteration"""
  dzprev::uType
  Δgzs::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  droptol::D
end
