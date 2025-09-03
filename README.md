# D layDiffEq.jl

[![Build Status](https://github.com/SciML/DelayDiffEq.jl/actions/workflows/Tests.yml/badge.svg?branch=master)](https://github.com/SciML/DelayDiffEq.jl/actions/workflows/Tests.yml?query=branch%3Amaster)
[![Coverage Status](https://coveralls.io/repos/SciML/DelayDiffEq.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/SciML/DelayDiffEq.jl?branch=master)
[![codecov](https://codecov.io/gh/SciML/DelayDiffEq.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/DelayDiffEq.jl)

DelayDiffEq.jl is a component package in the DifferentialEquations ecosystem. It holds the
delay differential equation solvers and utilities. It is built on top of OrdinaryDiffEq
to extend those solvers for delay differential equations. While completely independent
and usable on its own, users interested in using this
functionality should check out [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## API

DelayDiffEq.jl is part of the JuliaDiffEq common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes a DelayDiffEq.jl algorithm to `solve`. For example, we can solve the [DDE tutorial from the documentation](https://diffeq.sciml.ai/stable/tutorials/dde_example/) using the `MethodOfSteps(Tsit5())` algorithm:


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

For the list of available solvers, please refer to the [DifferentialEquations.jl DDE Solvers page](https://diffeq.sciml.ai/stable/solvers/dde_solve/). For options for the `solve` command, see the [common solver options page](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).

## Citing

If you use DelayDiffEq.jl in your work, please cite the following:

```bib
@article{DifferentialEquations.jl-2017,
 author = {Rackauckas, Christopher and Nie, Qing},
 doi = {10.5334/jors.151},
 journal = {The Journal of Open Research Software},
 keywords = {Applied Mathematics},
 note = {Exported from https://app.dimensions.ai on 2019/05/05},
 number = {1},
 pages = {},
 title = {DifferentialEquations.jl – A Performant and Feature-Rich Ecosystem for Solving Differential Equations in Julia},
 url = {https://app.dimensions.ai/details/publication/pub.1085583166 and http://openresearchsoftware.metajnl.com/articles/10.5334/jors.151/galley/245/download/},
 volume = {5},
 year = {2017}
}

@article{widmann2022delaydiffeq,
  title={DelayDiffEq: Generating Delay Differential Equation Solvers via Recursive Embedding of Ordinary Differential Equation Solvers},
  author={Widmann, David and Rackauckas, Chris},
  journal={arXiv preprint arXiv:2208.12879},
  year={2022}
}
```
