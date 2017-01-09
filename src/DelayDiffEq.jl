module DelayDiffEq

using DiffEqBase, OrdinaryDiffEq, DataStructures, RecursiveArrayTools
using Base.Test

import OrdinaryDiffEq: initialize!, perform_step!, loopfooter!,
       loopheader!, alg_order, handle_tstop!, ODEIntegrator, savevalues!

import DiffEqBase: solve, solve!, init


include("integrator_type.jl")
include("integrator_interface.jl")
include("history_function.jl")
include("algorithms.jl")
include("alg_utils.jl")
include("solve.jl")

export init

export initialize!, loopheader!, perform_step!, loopfooter!, handle_tstop!, postamble!

export HistoryFunction

export MethodOfSteps

end # module
