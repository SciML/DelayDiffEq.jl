using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, DataStructures
using Base.Test

import OrdinaryDiffEq: alg_order

const τ = 1
lags = [τ]
f = function (t,u,h)
  - h(t-τ)
end
h = (t) -> 0.0

function analytic(t,u0)
  t = t-1
  if t<0
    return u0
  elseif t<1
    return u0 - t*u0
  elseif t<2
    return 0.5*(3u0 -4*t*u0 + t^2 * u0)
  elseif t<3
    return 1/6*(17u0 -24t*u0 + 9t^2*u0 -t^3 * u0)
  elseif t<4
    return 1/24*(149u0 - 204*t*u0 + 90t^2*u0 - 16t^3*u0 + t^4*u0)
  elseif t<5
    return 1/120*(1769u0 - 2300t*u0 + 1090t^2*u0 - 240t^3*u0 + 25t^4*u0 - t^5*u0)
  elseif t<6
    return 1/720*(26239u0 - 32550t*u0 + 15915t^2*u0 - 3940t^3*u0 + 525t^4*u0 - 36t^5*u0 + t^6*u0)
  elseif t<7
    return (463609u0 - 554442t*u0 + 274701t^2*u0 - 72940t^3*u0 + 11235t^4*u0 - 1008t^5*u0 + 49t^6*u0 - t^7*u0)/5040
  elseif t<8
    return (9473673u0 - 11023880t*u0 + 5491780t^2*u0 - 1524712t^3*u0 + 257950t^4*u0 - 27272t^5*u0 + 1764t^6*u0 - 64t^7*u0 + t^8*u0)/40320
  elseif t<9
    return (219480785u0 - 250209864t*u0 + 124923492t^2*u0 - 35742504t^3*u0 +6450318t^4*u0 - 761544t^5*u0 + 58884t^6*u0 - 2880t^7*u0 + 81t^8*u0 - t^9*u0)/362880
  elseif t<=10
    return (219480785*u0 - 250209864*t*u0 + 124923492*t^2*u0 - 35742504t^3*u0 + 6450318t^4*u0 - 761544t^5*u0 + 58884t^6*u0 - 2880t^7*u0 + 81t^8*u0 - t^9*u0)/362880
  else
    error("This analytical solution is only valid on [-infty,11]")
  end
end


prob = DDETestProblem(f,h,1.0,analytic,lags,(0.0,10.0);iip=DiffEqBase.isinplace(f,4))
alg = MethodOfSteps(BS3();constrained=true)

dde_int = init(prob,alg;dt=0.1)
sol = solve!(dde_int)

h = (t) -> [0.0]
prob = DDEProblem(f,h,[1.0],lags,(0.0,10.0);iip=DiffEqBase.isinplace(f,4))
dde_int = init(prob,alg;dt=0.1)
sol = solve!(dde_int)

f = function (t,u,h,du)
  du[1] = - h(t-τ)[1]
end
h = (t) -> [0.0]
prob = DDEProblem(f,h,[1.0],lags,(0.0,10.0);iip=DiffEqBase.isinplace(f,4))
dde_int = init(prob,alg;dt=0.1)
sol = solve!(dde_int)


using Plots
plot(sol,ylim=[-0.51,1.05])

plott = linspace(0,10,100)
analytic_sol = [analytic(t,1.0) for t in linspace(0,10,100)]
plot!(plott,analytic_sol)

using RecursiveArrayTools
num_sol = vecvec_to_mat(sol(plott))
maximum(num_sol-analytic_sol)
