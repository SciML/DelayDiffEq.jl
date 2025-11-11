module DelayDiffEq

import Reexport: @reexport, Reexport
import OrdinaryDiffEqCore, OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqDifferentiation,
    OrdinaryDiffEqRosenbrock
import OrdinaryDiffEqDefault: OrdinaryDiffEqDefault
import OrdinaryDiffEqFunctionMap: OrdinaryDiffEqFunctionMap
@reexport using OrdinaryDiffEq

using DataStructures: BinaryMinHeap
using LinearAlgebra: opnorm, I
using Logging: @logmsg
using RecursiveArrayTools: copyat_or_push!, recursivecopy, recursivecopy!,
    recursive_bottom_eltype, recursive_unitless_bottom_eltype,
    recursive_unitless_eltype
using ForwardDiff: ForwardDiff

import ArrayInterface
import SimpleNonlinearSolve
import SymbolicIndexingInterface as SII

using SciMLBase: AbstractDDEAlgorithm, AbstractDDEIntegrator, AbstractODEIntegrator,
    DEIntegrator

using Base: deleteat!
import FastBroadcast: @..

using OrdinaryDiffEqNonlinearSolve: NLAnderson, NLFunctional
using OrdinaryDiffEqCore: AbstractNLSolverCache, SlowConvergence,
    alg_extrapolates, alg_maximum_order, initialize!, ODEVerbosity
using OrdinaryDiffEqRosenbrock: RosenbrockMutableCache
using OrdinaryDiffEqFunctionMap: FunctionMap
# using OrdinaryDiffEqDifferentiation: resize_grad_config!, resize_jac_config!

# Explicit imports for functions currently coming through @reexport using OrdinaryDiffEq
using OrdinaryDiffEqCore: AutoSwitch, CompositeAlgorithm
using OrdinaryDiffEq: OrdinaryDiffEq
using OrdinaryDiffEqDefault: DefaultODEAlgorithm
using SciMLBase: CallbackSet, DAEProblem, DDEProblem, DESolution, ODEProblem, ReturnCode,
    VectorContinuousCallback, addat!, addat_non_user_cache!,
    deleteat_non_user_cache!,
    du_cache, full_cache, get_tmp_cache, isinplace,
    reeval_internals_due_to_modification!,
    reinit!, resize_non_user_cache!, savevalues!, u_cache, user_cache,
    step!, terminate!, u_modified!, get_proposed_dt, set_proposed_dt!,
    auto_dt_reset!,
    add_tstop!, add_saveat!, get_du, get_du!, addsteps!,
    change_t_via_interpolation!, isadaptive
using DiffEqBase: initialize!
import DiffEqBase
using ConcreteStructs: @concrete
using SciMLLogging: AbstractVerbositySpecifier, AbstractVerbosityPreset, AbstractMessageLevel, None, Minimal, Standard, Detailed, All,
                    Silent, DebugLevel, InfoLevel, WarnLevel, ErrorLevel, @SciMLMessage

import SciMLBase

export Discontinuity, MethodOfSteps

include("discontinuity_type.jl")
include("functionwrapper.jl")
include("verbosity.jl")

include("integrators/type.jl")
include("integrators/utils.jl")
include("integrators/interface.jl")

include("fpsolve/type.jl")
include("fpsolve/fpsolve.jl")
include("fpsolve/utils.jl")
include("fpsolve/functional.jl")
include("cache_utils.jl")
include("interpolants.jl")
include("history_function.jl")
include("algorithms.jl")
include("track.jl")
include("alg_utils.jl")
include("solve.jl")
include("utils.jl")

# Default solver for DDEProblems
function SciMLBase.__solve(prob::DDEProblem; kwargs...)
    return SciMLBase.solve(prob, MethodOfSteps(DefaultODEAlgorithm()); kwargs...)
end
function SciMLBase.__init(prob::DDEProblem; kwargs...)
    return DiffEqBase.init(prob, MethodOfSteps(DefaultODEAlgorithm()); kwargs...)
end

end # module
