immutable HistoryFunction{F1,F2,F3} <: Function
  h::F1
  sol::F2
  integrator::F3
end

function (f::HistoryFunction)(t,deriv::Type=Val{0},idxs=nothing)
  @inbounds if t < f.sol.t[1]
    if typeof(idxs) <: Void
      return f.h(t)
    else
      return f.h(t,idxs)
    end
  elseif t <= f.sol.t[end] # Put equals back
    return f.sol.interp(t,idxs,deriv)
  else
    if typeof(idxs) <: Void
      return OrdinaryDiffEq.current_interpolant(t,f.integrator,eachindex(f.integrator.uprev),deriv)
    else
      return OrdinaryDiffEq.current_interpolant(t,f.integrator,idxs,deriv)
    end
  end
end

function (f::HistoryFunction)(val,t,deriv::Type=Val{0},idxs=nothing)
  @inbounds if t < f.sol.t[1]
    if typeof(idxs) <: Void
      return f.h(val,t)
    else
      return f.h(val,t,idxs)
    end
  elseif t <= f.sol.t[end] # Put equals back
    return f.sol.interp(val,t,idxs,deriv)
  else
    if typeof(idxs) <: Void
      return OrdinaryDiffEq.current_interpolant!(val,t,f.integrator,eachindex(f.integrator.uprev),deriv)
    else
      return OrdinaryDiffEq.current_interpolant!(val,t,f.integrator,idxs,deriv)
    end
  end
end
