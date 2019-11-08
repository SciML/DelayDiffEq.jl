# solver
mutable struct FPSolver{algType,iip,uTolType,C<:DiffEqBase.AbstractNLSolverCache} <: OrdinaryDiffEq.AbstractNLSolver{algType,iip}
  alg::algType
  κ::uTolType
  fast_convergence_cutoff::uTolType
  ηold::uTolType
  iter::Int
  maxiters::Int
  status::OrdinaryDiffEq.NLStatus
  cache::C
end

# caches
struct FPFunctionalCache{uType,uNoUnitsType} <: DiffEqBase.AbstractNLSolverCache
  atmp::uNoUnitsType
  dz::uType
end

struct FPFunctionalConstantCache <: DiffEqBase.AbstractNLSolverCache end

mutable struct FPAndersonCache{uType,uNoUnitsType,uEltypeNoUnits,D} <: DiffEqBase.AbstractNLSolverCache
  atmp::uNoUnitsType
  dz::uType
  dzold::uType
  z₊old::uType
  Δz₊s::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  aa_start::Int
  droptol::D
end

mutable struct FPAndersonConstantCache{uType,uEltypeNoUnits,D} <: DiffEqBase.AbstractNLSolverCache
  dz::uType
  dzold::uType
  z₊old::uType
  Δz₊s::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  aa_start::Int
  droptol::D
end
