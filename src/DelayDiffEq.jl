module DelayDiffEq

using Reexport
@reexport using OrdinaryDiffEq

using DataStructures
using LinearAlgebra
using Logging
using Printf
using RecursiveArrayTools
using SimpleUnPack

import ArrayInterface
import SimpleNonlinearSolve

using DiffEqBase: AbstractDDEAlgorithm, AbstractDDEIntegrator, AbstractODEIntegrator,
                  DEIntegrator, AbstractDDEProblem

using DiffEqBase: @..

if isdefined(OrdinaryDiffEq, :FastConvergence)
    using OrdinaryDiffEq: FastConvergence, Convergence, SlowConvergence,
                          VerySlowConvergence,
                          Divergence, AbstractNLSolverCache, NLNewton,
                          NLAnderson, NLFunctional
else
    using DiffEqBase: FastConvergence, Convergence, SlowConvergence, VerySlowConvergence,
                      Divergence, AbstractNLSolverCache
end
using OrdinaryDiffEq: RosenbrockMutableCache

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
