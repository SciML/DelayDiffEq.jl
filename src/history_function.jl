immutable HistoryFunction{F1,F2,F3} <: Function
  h::F1
  sol::F2
  integrator::F3
end

function (f::HistoryFunction)(t)
  if t < f.sol.t[1]
    return f.h(t)
  elseif t <= f.sol.t[end] # Put equals back
    return f.sol.interp(t)
  else
    return f.integrator(t)
  end
end
