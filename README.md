# DelayDiffEq.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/DelayDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DelayDiffEq.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ool429apgp28p71x?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/delaydiffeq-jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/DelayDiffEq.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/DelayDiffEq.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/DelayDiffEq.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/DelayDiffEq.jl?branch=master)

DelayDiffEq.jl is a component package in the DifferentialEquations ecosystem. It holds the
delay differential equation solvers and utilities. It is built on top of OrdinaryDiffEq
to extend those solvers for delay differential equations. While completely independent
and usable on its own, users interested in using this
functionality should check out [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

## API

DelayDiffEq.jl is part of the JuliaDiffEq common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes a DelayDiffEq.jl algorithm to `solve`. For example, we can solve the [DDE tutorial from the documentation](http://docs.juliadiffeq.org/latest/tutorials/dde_example.html) using the `MethodOfSteps(Tsit5())` algorithm:

```julia
using DelayDiffEq
const p0 = 0.2; const q0 = 0.3; const v0 = 1; const d0 = 5
const p1 = 0.2; const q1 = 0.3; const v1 = 1; const d1 = 1
const d2 = 1; const beta0 = 1; const beta1 = 1; const tau = 1
function bc_model(du,u,h,p,t)
  du[1] = (v0/(1+beta0*(h(p, t-tau)[3]^2))) * (p0 - q0)*u[1] - d0*u[1]
  du[2] = (v0/(1+beta0*(h(p, t-tau)[3]^2))) * (1 - p0 + q0)*u[1] +
          (v1/(1+beta1*(h(p, t-tau)[3]^2))) * (p1 - q1)*u[2] - d1*u[2]
  du[3] = (v1/(1+beta1*(h(p, t-tau)[3]^2))) * (1 - p1 + q1)*u[2] - d2*u[3]
end
lags = [tau]
h(p, t) = ones(3)
tspan = (0.0,10.0)
u0 = [1.0,1.0,1.0]
prob = DDEProblem(bc_model,u0,h,tspan,constant_lags = lags)
alg = MethodOfSteps(Tsit5())
sol = solve(prob,alg)
using Plots; plot(sol)
```

Both constant and state-dependent lags are supported. Interfacing with OrdinaryDiffEq.jl for implicit methods for stiff equations is also supported.

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl DDE Solvers page](http://docs.juliadiffeq.org/latest/solvers/dde_solve.html). For options for the `solve` command, see the [common solver options page](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html).
