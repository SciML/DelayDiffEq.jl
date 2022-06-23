## SciMLBase Trait Definitions

function SciMLBase.isautodifferentiable(alg::AbstractMethodOfStepsAlgorithm)
    SciMLBase.isautodifferentiable(alg.alg)
end
function SciMLBase.allows_arbitrary_number_types(alg::AbstractMethodOfStepsAlgorithm)
    SciMLBase.allows_arbitrary_number_types(alg.alg)
end
function SciMLBase.allowscomplex(alg::AbstractMethodOfStepsAlgorithm)
    SciMLBase.allowscomplex(alg.alg)
end
SciMLBase.isdiscrete(alg::AbstractMethodOfStepsAlgorithm) = SciMLBase.isdiscrete(alg.alg)
SciMLBase.isadaptive(alg::AbstractMethodOfStepsAlgorithm) = SciMLBase.isadaptive(alg.alg)

## DelayDiffEq Internal Traits

function isconstrained(alg::AbstractMethodOfStepsAlgorithm{constrained}) where {constrained}
    constrained
end
OrdinaryDiffEq.uses_uprev(alg::AbstractMethodOfStepsAlgorithm, adaptive) = true

function OrdinaryDiffEq.isimplicit(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.isimplicit(alg.alg)
end
function OrdinaryDiffEq.isdtchangeable(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.isdtchangeable(alg.alg)
end
function OrdinaryDiffEq.ismultistep(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.ismultistep(alg.alg)
end
function OrdinaryDiffEq.isautoswitch(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.isautoswitch(alg.alg)
end
function OrdinaryDiffEq.get_chunksize(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.get_chunksize(alg.alg)
end
function OrdinaryDiffEq.get_chunksize_int(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.get_chunksize_int(alg.alg)
end
function OrdinaryDiffEq.alg_autodiff(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.alg_autodiff(alg.alg)
end
function OrdinaryDiffEq.alg_difftype(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.alg_difftype(alg.alg)
end
function OrdinaryDiffEq.standardtag(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.standardtag(alg.alg)
end
function OrdinaryDiffEq.concrete_jac(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.concrete_jac(alg.alg)
end
function OrdinaryDiffEq.alg_extrapolates(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.alg_extrapolates(alg.alg)
end
function OrdinaryDiffEq.alg_order(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.alg_order(alg.alg)
end
function OrdinaryDiffEq.alg_maximum_order(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.alg_maximum_order(alg.alg)
end
function OrdinaryDiffEq.alg_adaptive_order(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEq.alg_adaptive_order(alg.alg)
end

"""
    iscomposite(alg)

Return if algorithm `alg` is a composite algorithm.
"""
iscomposite(alg) = false
iscomposite(::OrdinaryDiffEq.OrdinaryDiffEqCompositeAlgorithm) = true
iscomposite(alg::AbstractMethodOfStepsAlgorithm) = iscomposite(alg.alg)
