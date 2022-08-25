using CSV
using DataFrames
using DelayDiffEq

@doc raw"""
    mackey_glass()

Solve the Mackey-Glass equation on the time interval `[0, 600]`.

The Mackey-Glass equation

```math
\begin{aligned}
u'(t) &= 2 u(t - 2) / (1 + u(t - 2)^9.65) - u(t), \qquad t \geq 0,\\
u(t) &= 0.5, \qquad t \leq 0.
\end{aligned}
```

is a delay-differential equation model of circulating blood cells that exhibits a chaotic solution.

## Reference

Glass, L., and M. C. Mackey. (1979). Pathological conditions resulting from instabilities in physiological control systems, Annals of the New York Academy of Sciences, 316(1), pp. 214-235.
"""
function mackey_glass()
    # Define DDE problem
    function mackey_glass_equation(u, h, p, t)
        utau = h(p, t - 2)
        return 2 * utau / (1 + utau^9.65) - u
    end
    mackey_glass_history(p, t) = 0.5
    tspan = (0.0, 600.0)
    prob = DDEProblem(mackey_glass_equation,
                      mackey_glass_history, tspan;
                      constant_lags = (2,))

    # Compute solution
    alg = MethodOfSteps(Rodas5())
    sol = solve(prob, alg; reltol = 1e-14,
                abstol = 1e-14)

    return sol
end

"""
    mackey_glass(io)

Write the solution computed by [`waltman()`](@ref) for time points `300:0.1:600` together with the delayed state to `io` in the CSV file format, and return the solution object.

The CSV file consists of columns `t` (time), `x` (state), and `xtau` (delayed state).

!!! note
    If `io` is a file name, its directory is created if it does not exist yet.
"""
function mackey_glass(io)
    # Solve the DDE problem
    sol = mackey_glass()

    # Create dataframe of the solution at desired time points
    ts = 300:0.1:600
    df = DataFrame(; t = ts, x = sol(ts).u,
                   xtau = sol(ts .- 2).u)

    # Write data frame in CSV file format
    io isa String && mkpath(dirname(io))
    CSV.write(io, df)

    return sol
end
