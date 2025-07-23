# DelayDiffEq.jl: Delay Differential Equation Solvers

DelayDiffEq.jl is a component package of the DifferentialEquations.jl ecosystem for solving delay differential equations (DDEs). It provides the core algorithms and method of steps implementations for solving DDEs with both constant and state-dependent delays.

## Features

- Method of steps algorithms with automatic step size control
- Support for both constant and state-dependent delays
- Compatible with all ODE solvers from OrdinaryDiffEq.jl
- Automatic discontinuity tracking for accurate solutions
- Full compatibility with the DifferentialEquations.jl common interface

## Installation

To install DelayDiffEq.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("DelayDiffEq")
```

## Quick Example

```julia
using DelayDiffEq, Plots

# Define a DDE with constant delay
function bc_model(du,u,h,p,t)
    du[1] = 1.1/(1 + sqrt(10)*(h(p, t-20)[1])^(5/4)) - 10*u[1]/(1 + 40*u[2])
    du[2] = 100*u[1]/(1 + 40*u[2]) - 2.43*u[2]
end

h(p, t) = ones(2)
tspan = (0.0,100.0)
u0 = [1.05767027/3, 1.030713491/3]

prob = DDEProblem(bc_model,u0,h,tspan; constant_lags=[20.0])
sol = solve(prob,MethodOfSteps(Tsit5()))
plot(sol)
```

## Getting Started

For more examples and tutorials, see the [DifferentialEquations.jl DDE tutorial](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/dde_example/).

For the list of available algorithms and their properties, see the [Solver API](@ref) documentation.