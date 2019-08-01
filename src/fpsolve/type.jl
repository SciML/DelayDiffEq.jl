abstract type AbstractFPSolverAlgorithm end
abstract type AbstractFPSolverCache end

@enum FPStatus::Int8 begin
  FastConvergence     = 2
  Convergence         = 1
  SlowConvergence     = 0
  VerySlowConvergence = -1
  Divergence          = -2
end

# solver
mutable struct FPSolver{C<:AbstractFPSolverCache,uTolType,K,C1}
  ηold::uTolType
  κ::K
  max_iter::Int
  fp_iters::Int
  status::FPStatus
  fast_convergence_cutoff::C1
  cache::C
end

struct FPNone end

# algorithms
struct FPFunctional{K,C} <: AbstractFPSolverAlgorithm
  κ::K
  fast_convergence_cutoff::C
  max_iter::Int
end

FPFunctional(; κ=1//100, max_iter=10, fast_convergence_cutoff=1//5) =
  FPFunctional(κ, fast_convergence_cutoff, max_iter)

# caches
struct FPFunctionalCache{uType} <: AbstractFPSolverCache
  resid::uType
end

struct FPFunctionalConstantCache <: AbstractFPSolverCache end
