isconstrained(alg::AbstractMethodOfStepsAlgorithm{constrained}) where constrained =
    constrained

OrdinaryDiffEq.alg_maximum_order(alg::AbstractMethodOfStepsAlgorithm) =
  OrdinaryDiffEq.alg_maximum_order(alg.alg)
OrdinaryDiffEq.alg_extrapolates(alg::AbstractMethodOfStepsAlgorithm) =
  OrdinaryDiffEq.alg_extrapolates(alg.alg)

OrdinaryDiffEq.uses_uprev(alg::AbstractMethodOfStepsAlgorithm, adaptive) = true

DiffEqBase.isadaptive(alg::AbstractMethodOfStepsAlgorithm) =
  DiffEqBase.isadaptive(alg.alg)

"""
    iscomposite(alg)

Return if algorithm `alg` is a composite algorithm.
"""
iscomposite(alg) = false
iscomposite(::OrdinaryDiffEq.OrdinaryDiffEqCompositeAlgorithm) = true
iscomposite(alg::AbstractMethodOfStepsAlgorithm) = iscomposite(alg.alg)

