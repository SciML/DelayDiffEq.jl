using DelayDiffEq
using SciMLSensitivity, ForwardDiff, Zygote

# Lotka-Volterra model with delay
function lv!(du, u, h, p, t)
    x, y = u
    α, β, δ, γ = p
    xτ = h(p, t - 0.1)[1]
    du[1] = dx = (α - β * y) * xτ
    du[2] = dy = (δ * x - γ) * y
    return nothing
end

# Constant history function
history_lv(p, t) = ones(2)

"""
    predict_lv(p)

Compute the solution for the Lotka-Volterra model with delay, defined by [`lv!`](@ref)
and `history_lv`, at time points `t ∈ 0:0.1:5` for parameters `p = [α, β, δ, γ]` as
an array of dimension state x time.
"""
function predict_lv(p)
    # Define the DDE problem
    prob = DDEProblem(lv!, history_lv, (0.0, 5.0),
                      p;
                      constant_lags = (0.1,))

    # Solve the problem and save the solution at
    # time points 0, 0.1, 0.2, ..., 5
    sol = solve(prob, MethodOfSteps(Tsit5());
                saveat = 0.1)

    return Array(sol)
end

"""
    jacobian_lv()

Compute the Jacobian of [`predict_lv`](@ref) for parameters `p = [2.2, 1.0, 2.0, 0.4]` with
ForwardDiff (forward-mode) and Zygote (reverse-mode) AD backends.
"""
function jacobian_lv()
    # Parameters
    p = [2.2, 1.0, 2.0, 0.4]

    # Compute Jacobian
    jacFD = ForwardDiff.jacobian(predict_lv, p)
    jacZygote = Zygote.jacobian(predict_lv, p)[1]

    return (; forwarddiff = jacFD,
            zygote = jacZygote)
end

## Separate documentation
## (fixes issues with listing package in LaTeX)

@doc raw"""
    lv!(du, u, h, p, t)

Compute the derivatives of a Lotka-Volterra model with delay at time `t` and state 
`u = [x, y]` with parameters `p = [α, β, δ, γ]` and history function `h`, and write them to
`du = [dx, dy]`.

The delay differential equation model is given by
```math
\begin{aligned}
x'(t) &= (α - β y(t)) x(t - 0.1),\\
y'(t) &= (δ x(t) - γ) y(t).
\end{aligned}
```
""" lv!
