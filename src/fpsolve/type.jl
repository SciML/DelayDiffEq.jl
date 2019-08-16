# solver
mutable struct FPSolver{algType<:DiffEqBase.AbstractNLSolverAlgorithm,IIP,uTolType,C}
  alg::algType
  ηold::uTolType
  fp_iters::Int
  status::DiffEqBase.NLStatus
  cache::C
end

# caches
struct FPFunctionalCache{uType,uNoUnitsType} <: DiffEqBase.AbstractNLSolverCache
  du::uType
  atmp::uNoUnitsType
end

struct FPFunctionalConstantCache <: DiffEqBase.AbstractNLSolverCache end

struct FPAndersonCache{uType,uNoUnitsType,gsType,QType,RType,gType} <: DiffEqBase.AbstractNLSolverCache
  du::uType
  duold::uType
  uold::uType
  atmp::uNoUnitsType
  Δus::gsType
  Q::QType
  R::RType
  γs::gType
end

struct FPAndersonConstantCache{gsType,QType,RType,gType} <: DiffEqBase.AbstractNLSolverCache
  Δus::gsType
  Q::QType
  R::RType
  γs::gType
end
