if haskey(ENV,"GROUP")
   group = ENV["GROUP"]
else
   group = "All"
end

is_APPVEYOR = ( Sys.iswindows() && haskey(ENV,"APPVEYOR") )

if group == "All" || group == "Interface"
  tests = ["discontinuities.jl",
           "history_function.jl",
           "parameters.jl",
           "retcode.jl",
           "composite_solution.jl",
           "dependent_delays.jl",
           "reinit.jl",
           "saveat.jl",
           "save_idxs.jl",
           "events.jl",
           "units.jl",
           "unique_times.jl"
           ]
   for test in tests
       @time include(test)
   end
elseif group == "All" || group == "Integrators"
  tests =   ["constrained.jl",
             "unconstrained.jl",
             "residual_control.jl",
             "lazy_interpolants.jl",
             "sdirk_integrators.jl",
             "rosenbrock_integrators.jl"
             ]
   for test in tests
       @time include(test)
   end
end
