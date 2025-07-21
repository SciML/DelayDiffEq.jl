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
OrdinaryDiffEqCore.uses_uprev(alg::AbstractMethodOfStepsAlgorithm, adaptive) = true

function OrdinaryDiffEqCore.isimplicit(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.isimplicit(alg.alg)
end
function OrdinaryDiffEqCore.isdtchangeable(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.isdtchangeable(alg.alg)
end
function OrdinaryDiffEqCore.ismultistep(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.ismultistep(alg.alg)
end
function OrdinaryDiffEqCore.isautoswitch(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.isautoswitch(alg.alg)
end
function OrdinaryDiffEqCore.get_chunksize(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.get_chunksize(alg.alg)
end
function OrdinaryDiffEqCore.get_chunksize_int(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.get_chunksize_int(alg.alg)
end
function OrdinaryDiffEqCore.alg_autodiff(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.alg_autodiff(alg.alg)
end
function OrdinaryDiffEqCore.alg_difftype(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.alg_difftype(alg.alg)
end
function OrdinaryDiffEqCore.standardtag(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.standardtag(alg.alg)
end
function OrdinaryDiffEqCore.concrete_jac(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.concrete_jac(alg.alg)
end
function OrdinaryDiffEqCore.alg_extrapolates(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.alg_extrapolates(alg.alg)
end
function SciMLBase.alg_order(alg::AbstractMethodOfStepsAlgorithm)
    SciMLBase.alg_order(alg.alg)
end
function OrdinaryDiffEqCore.alg_maximum_order(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.alg_maximum_order(alg.alg)
end
function OrdinaryDiffEqCore.alg_adaptive_order(alg::AbstractMethodOfStepsAlgorithm)
    OrdinaryDiffEqCore.alg_adaptive_order(alg.alg)
end

"""
    iscomposite(alg)

Return if algorithm `alg` is a composite algorithm.
"""
iscomposite(alg) = false
iscomposite(::OrdinaryDiffEqCore.OrdinaryDiffEqCompositeAlgorithm) = true
iscomposite(alg::AbstractMethodOfStepsAlgorithm) = iscomposite(alg.alg)

function DiffEqBase.prepare_alg(alg::MethodOfSteps, u0, p, prob)
    MethodOfSteps(DiffEqBase.prepare_alg(alg.alg, u0, p, prob);
        constrained = isconstrained(alg),
        fpsolve = alg.fpsolve)
end
