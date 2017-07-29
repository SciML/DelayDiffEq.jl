__precompile__()

module DelayDiffEq

using Reexport
@reexport using DiffEqBase

using OrdinaryDiffEq, DataStructures, RecursiveArrayTools, Combinatorics

import OrdinaryDiffEq: initialize!, perform_step!, loopfooter!, loopheader!, alg_order,
                       handle_tstop!, ODEIntegrator, savevalues!, handle_callback_modifiers!

import DiffEqBase: solve, solve!, init, resize!, u_cache, user_cache, du_cache, full_cache,
                   deleteat!, terminate!, u_modified!, get_proposed_dt, set_proposed_dt!

import OrdinaryDiffEq: Rosenbrock23Cache, Rosenbrock32Cache, ImplicitEulerCache,
                       TrapezoidCache

include("integrator_type.jl")
include("integrator_utils.jl")
include("integrator_interface.jl")
include("history_function.jl")
include("algorithms.jl")
include("callbacks.jl")
include("alg_utils.jl")
include("solve.jl")
include("utils.jl")

export init

export initialize!, loopheader!, perform_step!, loopfooter!, handle_tstop!, postamble!

export HistoryFunction

export MethodOfSteps

end # module
