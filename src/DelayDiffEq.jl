module DelayDiffEq

using DataStructures
using LinearAlgebra
using Logging
using Printf
using RecursiveArrayTools
using SimpleUnPack

import ArrayInterface
import SimpleNonlinearSolve
import SymbolicIndexingInterface as SII

using DiffEqBase: DiffEqBase, AbstractDDEAlgorithm, AbstractDDEIntegrator, AbstractODEIntegrator,
                  DEIntegrator, AbstractDDEProblem

using FastBroadcast
using Reexport
@reexport using DiffEqBase

import OrdinaryDiffEqCore, OrdinaryDiffEqNonlinearSolve, OrdinaryDiffEqDifferentiation, OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore: CompositeAlgorithm, AutoSwitch
using OrdinaryDiffEqFunctionMap: FunctionMap
using OrdinaryDiffEqRosenbrock: RosenbrockMutableCache
using OrdinaryDiffEqNonlinearSolve: NLNewton, NLAnderson, NLFunctional, AbstractNLSolverCache,
    FastConvergence, Convergence, SlowConvergence, VerySlowConvergence, Divergence

import SciMLBase

export DDEProblem

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
