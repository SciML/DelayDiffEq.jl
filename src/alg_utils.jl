isconstrained(alg::AbstractMethodOfStepsAlgorithm{constrained}) where constrained =
    constrained
OrdinaryDiffEq.uses_uprev(alg::AbstractMethodOfStepsAlgorithm, adaptive) = true

OrdinaryDiffEq.isimplicit(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.isimplicit(alg.alg)
OrdinaryDiffEq.isdtchangeable(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.isdtchangeable(alg.alg)
OrdinaryDiffEq.ismultistep(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.ismultistep(alg.alg)
OrdinaryDiffEq.isadaptive(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.isadaptive(alg.alg)
OrdinaryDiffEq.isautoswitch(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.isautoswitch(alg.alg)
OrdinaryDiffEq.get_chunksize(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.get_chunksize(alg.alg)
OrdinaryDiffEq.get_chunksize_int(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.get_chunksize_int(alg.alg)
OrdinaryDiffEq.alg_autodiff(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.alg_autodiff(alg.alg)
OrdinaryDiffEq.alg_difftype(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.alg_difftype(alg.alg)
OrdinaryDiffEq.standardtag(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.standardtag(alg.alg)
OrdinaryDiffEq.concrete_jac(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.concrete_jac(alg.alg)
OrdinaryDiffEq.alg_extrapolates(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.alg_extrapolates(alg.alg)
OrdinaryDiffEq.alg_order(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.alg_order(alg.alg)
OrdinaryDiffEq.alg_maximum_order(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.alg_maximum_order(alg.alg)
OrdinaryDiffEq.alg_adaptive_order(alg::AbstractMethodOfStepsAlgorithm) = OrdinaryDiffEq.alg_adaptive_order(alg.alg)

"""
    iscomposite(alg)

Return if algorithm `alg` is a composite algorithm.
"""
iscomposite(alg) = false
iscomposite(::OrdinaryDiffEq.OrdinaryDiffEqCompositeAlgorithm) = true
iscomposite(alg::AbstractMethodOfStepsAlgorithm) = iscomposite(alg.alg)
