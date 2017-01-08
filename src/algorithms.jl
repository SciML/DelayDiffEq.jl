abstract DelayDiffEqAlgorithm <: AbstractDDEAlgorithm
abstract AbstractMethodOfStepsAlgorithm{constrained} <: DelayDiffEqAlgorithm

immutable MethodOfSteps{algType,constrained} <: AbstractMethodOfStepsAlgorithm{constrained}
  alg::algType
end

Base.@pure MethodOfSteps(alg;constrained=false) =  MethodOfSteps{typeof(alg),constrained}(alg)
