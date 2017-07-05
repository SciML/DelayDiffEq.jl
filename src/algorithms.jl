@compat abstract type DelayDiffEqAlgorithm <: AbstractDDEAlgorithm end
@compat abstract type AbstractMethodOfStepsAlgorithm{constrained} <: DelayDiffEqAlgorithm end

immutable MethodOfSteps{algType,AType,RType,NType,constrained} <: AbstractMethodOfStepsAlgorithm{constrained}
  alg::algType
  fixedpoint_abstol::AType
  fixedpoint_reltol::RType
  max_fixedpoint_iters::Int
  picardnorm::NType # anderson acceleration with nlsolve always uses infinite norm of residuals
  m::Int # controls history size of anderson acceleration (m=0 corresponds to simple fixed-point iteration)
end

Base.@pure MethodOfSteps(alg;constrained=false,
                             fixedpoint_abstol = nothing,
                             fixedpoint_reltol = nothing,
                             max_fixedpoint_iters = 10,
                             picardnorm = nothing,
                             m = 0) =
                             MethodOfSteps{typeof(alg),typeof(fixedpoint_abstol),typeof(fixedpoint_reltol),typeof(picardnorm),constrained}(alg,fixedpoint_abstol,fixedpoint_reltol,max_fixedpoint_iters,picardnorm,m)
