```@meta
CollapsedDocStrings = true
```

# Solver API

DelayDiffEq.jl provides delay differential equation solvers through the method of steps approach. The core algorithm wraps ODE solvers from OrdinaryDiffEq.jl to handle the delays.

## Method of Steps Algorithm

The primary algorithm for solving delay differential equations in DelayDiffEq.jl is `MethodOfSteps`, which implements the method of steps approach for solving DDEs.

```@docs
MethodOfSteps
```

### Algorithm Properties

- **Adaptive**: Inherits adaptivity from the underlying ODE solver
- **Order**: Depends on the chosen ODE algorithm
- **Dense Output**: Available when the underlying ODE solver supports it
- **State-Dependent Delays**: Supported through fixed-point iteration

## Usage

### Basic Usage

```julia
using DelayDiffEq, OrdinaryDiffEq

# Use with any ODE solver from OrdinaryDiffEq.jl
alg = MethodOfSteps(Tsit5())
sol = solve(prob, alg)
```

### Constrained vs Unconstrained

By default, `MethodOfSteps` is unconstrained, meaning it can take steps larger than the minimal delay using fixed-point iteration:

```julia
# Unconstrained (default) - can take larger steps
alg = MethodOfSteps(Tsit5())

# Constrained - steps limited to minimal delay
alg = MethodOfSteps(Tsit5(); constrained = true)
```

### Custom Fixed-Point Solver

For unconstrained problems, you can specify the fixed-point iteration method:

```julia
# Use custom fixed-point solver
alg = MethodOfSteps(Tsit5(); fpsolve = NLFunctional(; max_iter = 100))
```

## Recommended Algorithms

The choice of underlying ODE algorithm depends on your problem characteristics:

### Non-Stiff Problems

For non-stiff delay differential equations:

- **`MethodOfSteps(Tsit5())`**: Good general purpose solver (5th order)
- **`MethodOfSteps(BS3())`**: For lower accuracy requirements (3rd order)
- **`MethodOfSteps(Vern6())`**: For high accuracy requirements (6th order)
- **`MethodOfSteps(Vern9())`**: For very high accuracy requirements (9th order)

### Stiff Problems

For stiff delay differential equations:

- **`MethodOfSteps(Rosenbrock23())`**: Good for mildly stiff problems
- **`MethodOfSteps(Rodas4())`**: General purpose stiff solver
- **`MethodOfSteps(Rodas5())`**: High accuracy stiff solver
- **`MethodOfSteps(KenCarp4())`**: Good stability properties

### Low Storage Requirements

For problems with memory constraints:

- **`MethodOfSteps(SSPRK104())`**: Strong stability preserving, low storage
- **`MethodOfSteps(OwrenZen3())`**: 3rd order low storage
- **`MethodOfSteps(OwrenZen5())`**: 5th order low storage

## Algorithm Selection Guide

### Step 1: Determine Stiffness

First, determine if your DDE is stiff. Signs of stiffness include:
- Explicit methods require very small time steps
- The solution has multiple time scales
- There are rapid transients followed by slow dynamics

### Step 2: Choose Tolerance

- **High tolerance (>1e-2)**: Use lower order methods
- **Medium tolerance (1e-8 to 1e-2)**: Use standard 4-5 order methods  
- **Low tolerance (<1e-8)**: Use high order methods

### Step 3: Consider Problem Structure

- **Discontinuous forcing**: Use low order methods or constrained stepping
- **State-dependent delays**: May benefit from constrained stepping
- **Conservation properties needed**: Consider symplectic integrators

### Example Selection

```julia
# Non-stiff problem with medium accuracy
alg = MethodOfSteps(Tsit5())

# Stiff problem with medium accuracy
alg = MethodOfSteps(Rodas4())

# High accuracy non-stiff problem
alg = MethodOfSteps(Vern9())

# Problem with frequent discontinuities
alg = MethodOfSteps(BS3(); constrained = true)
```

## Performance Tips

1. **Use constrained stepping** when delays are much smaller than the timescale of the solution
2. **Choose higher order methods** for smooth problems with stringent accuracy requirements
3. **Use lower order methods** when the solution has many discontinuities
4. **Monitor the number of fixed-point iterations** for unconstrained problems with state-dependent delays

## See Also

- [OrdinaryDiffEq.jl Solver Documentation](https://docs.sciml.ai/OrdinaryDiffEq/stable/) for details on the underlying ODE solvers
- [DifferentialEquations.jl DDE Tutorial](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/dde_example/) for comprehensive examples