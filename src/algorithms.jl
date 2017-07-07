abstract type DelayDiffEqAlgorithm <: AbstractDDEAlgorithm end
abstract type AbstractMethodOfStepsAlgorithm{constrained} <: DelayDiffEqAlgorithm end

struct MethodOfSteps{algType,AType,RType,NType,constrained} <:
    AbstractMethodOfStepsAlgorithm{constrained}

    alg::algType
    fixedpoint_abstol::AType
    fixedpoint_reltol::RType
    fixedpoint_norm::NType
    max_fixedpoint_iters::Int
end

"""
    MethodOfSteps(alg; constrained::Bool=false, fixedpoint_abstol=nothing,
                  fixedpoint_reltol=nothing, fixedpoint_norm=nothing,
                  max_fixedpoint_iters::Int=10)

Construct an algorithm that solves delay differential equations by the method of steps,
where `alg` is an ODE algorithm from OrdinaryDiffEq.jl without lazy interpolation upon
which the calculation of steps is based.

If the algorithm is `constrained` only steps of size at most the minimal delay will be
taken. If it is unconstrained, fixed-point iteration is applied for step sizes that exceed
the minimal delay.

The absolute and relative tolerance of the fixed-point iterations can be set by
`fixedpoint_abstol` and `fixedpoint_reltol`, respectively, either as scalars or vectors.
Based on these tolerances error estimates are calculated during the fixed-point iterations
with a norm that may be specified as `fixedpoint_norm`. Fixed-point iterations are stopped
if the error estimate is less than 1 or after the maximal number `max_fixedpoint_iters` of
iteration steps.
"""
Base.@pure function MethodOfSteps(alg; constrained::Bool=false, fixedpoint_abstol=nothing,
                                  fixedpoint_reltol=nothing, fixedpoint_norm=nothing,
                                  max_fixedpoint_iters::Int=10)
    MethodOfSteps{typeof(alg),typeof(fixedpoint_abstol),typeof(fixedpoint_reltol),
                  typeof(fixedpoint_norm),constrained}(alg,fixedpoint_abstol,
                                                       fixedpoint_reltol,
                                                       fixedpoint_norm, max_fixedpoint_iters)
end
