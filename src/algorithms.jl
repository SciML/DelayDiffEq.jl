abstract type DelayDiffEqAlgorithm <: AbstractDDEAlgorithm end
abstract type AbstractMethodOfStepsAlgorithm{constrained} <: DelayDiffEqAlgorithm end

struct MethodOfSteps{algType,F,constrained} <: AbstractMethodOfStepsAlgorithm{constrained}
  alg::algType
  fpsolve::F
end

"""
    MethodOfSteps(alg; constrained = false, fpsolve = NLFunctional())

Construct an algorithm that solves delay differential equations by the method of steps,
where `alg` is an ODE algorithm from OrdinaryDiffEq.jl upon which the calculation of
steps is based.

If the algorithm is `constrained` only steps of size at most the minimal delay will be
taken. If it is unconstrained, fixed-point iteration `fpsolve` is applied for step sizes
that exceed the minimal delay.
"""
MethodOfSteps(alg; constrained = false, fpsolve = NLFunctional()) =
  MethodOfSteps{typeof(alg),typeof(fpsolve),constrained}(alg, fpsolve)
