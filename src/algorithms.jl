@compat abstract type DelayDiffEqAlgorithm <: AbstractDDEAlgorithm end
@compat abstract type AbstractMethodOfStepsAlgorithm{constrained} <: DelayDiffEqAlgorithm end

immutable MethodOfSteps{algType,AType,RType,NType,constrained} <: AbstractMethodOfStepsAlgorithm{constrained}
  alg::algType
  picardabstol::AType
  picardreltol::RType
  picardnorm::NType
  max_picard_iters::Int
end

Base.@pure MethodOfSteps(alg;constrained=false,
                             picardabstol = nothing,
                             picardreltol = nothing,
                             picardnorm   = nothing,
                             max_picard_iters = 10) =
                             MethodOfSteps{typeof(alg),typeof(picardabstol),typeof(picardreltol),typeof(picardnorm),constrained}(alg,picardabstol,picardreltol,picardnorm,max_picard_iters)
