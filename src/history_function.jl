immutable HistoryFunction{F1,F2,F3} <: Function
  h::F1
  sol::F2
  integrator::F3
end

function (f::HistoryFunction)(t,deriv::Type=Val{0};idxs=nothing)
  if t < f.sol.t[1]
    if typeof(idxs) <: Void
      return f.h(t)
    else
      return f.h(t,idxs)
    end
  elseif t <= f.sol.t[end] # Put equals back
    return f.sol.interp(t,idxs,deriv)
  else
    return OrdinaryDiffEq.current_interpolant(t,f.integrator,idxs,deriv)
  end
end

function (f::HistoryFunction)(val,t,deriv::Type=Val{0};idxs=nothing)
  if typeof(idxs) <: Void
    if t < f.sol.t[1]
      return f.h(val,t)
    elseif t <= f.sol.t[end] # Put equals back
      return f.sol.interp(val,t,idxs,deriv)
    else
      return f.integrator(val,t,deriv)
    end
  else
    if t < f.sol.t[1]
      return f.h(val,t)
    elseif t <= f.sol.t[end] # Put equals back
      return f.sol.interp(val,t,idxs,deriv)
    else
      return f.integrator(val,t,deriv;idxs=idxs)
    end
  end
end
