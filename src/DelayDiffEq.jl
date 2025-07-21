module DelayDiffEq

using Reexport: @reexport
import OrdinaryDiffEqCore, OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqDifferentiation, OrdinaryDiffEqRosenbrock
@reexport using OrdinaryDiffEq

using DataStructures: BinaryMinHeap
using LinearAlgebra: opnorm, I
using Logging: @logmsg
using RecursiveArrayTools: copyat_or_push!, recursivecopy, recursivecopy!, recursive_bottom_eltype, recursive_unitless_bottom_eltype, recursive_unitless_eltype
using SimpleUnPack: @unpack

import ArrayInterface
import SimpleNonlinearSolve
import SymbolicIndexingInterface as SII

using SciMLBase: AbstractDDEAlgorithm, AbstractDDEIntegrator, AbstractODEIntegrator,
                  DEIntegrator

using Base: deleteat!
import FastBroadcast: @..

using OrdinaryDiffEqNonlinearSolve: NLAnderson, NLFunctional, NLStatus, nlsolvefail, initial_η, apply_step!
using OrdinaryDiffEqCore: AbstractNLSolverCache, SlowConvergence, OrdinaryDiffEqCompositeAlgorithm, InterpolationData, 
                          DEOptions, AutoSwitchCache, alg_cache, get_fsalfirstlast, ode_determine_initdt,
                          _ode_addsteps!, default_controller, gamma_default, qmin_default, qmax_default,
                          qsteady_min_default, qsteady_max_default, initialize_tstops, initialize_saveat,
                          initialize_d_discontinuities, alg_extrapolates, alg_maximum_order, get_differential_vars,
                          _initialize_dae!, loopheader!, loopfooter!, perform_step!, initialize!, handle_dt!,
                          handle_tstop!, initialize_callbacks!, nlsolve_f, current_interpolant, current_interpolant!
using OrdinaryDiffEqRosenbrock: RosenbrockMutableCache
using OrdinaryDiffEqFunctionMap: FunctionMap
# using OrdinaryDiffEqDifferentiation: resize_grad_config!, resize_jac_config!

# Explicit imports for functions currently coming through @reexport using OrdinaryDiffEq
using OrdinaryDiffEqCore: AutoSwitch, CompositeAlgorithm
using OrdinaryDiffEq: OrdinaryDiffEq
using SciMLBase: CallbackSet, DAEProblem, DDEProblem, DESolution, ODEProblem, ReturnCode,
                 VectorContinuousCallback, addat!, addat_non_user_cache!, deleteat_non_user_cache!,
                 du_cache, full_cache, get_tmp_cache, isinplace, reeval_internals_due_to_modification!,
                 reinit!, resize_non_user_cache!, savevalues!, u_cache, user_cache,
                 step!, terminate!, u_modified!, get_proposed_dt, set_proposed_dt!, auto_dt_reset!,
                 add_tstop!, add_saveat!, get_du, get_du!, addsteps!,
                 change_t_via_interpolation!, isadaptive
using DiffEqBase: initialize!
import DiffEqBase

import SciMLBase

export Discontinuity, MethodOfSteps

include("discontinuity_type.jl")
include("functionwrapper.jl")

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

end # module
