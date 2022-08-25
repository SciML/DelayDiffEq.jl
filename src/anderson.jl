using DelayDiffEq
using DDEProblemLibrary

@doc raw"""
    compare_anderson()

Solve the delay-differential equation problem

```math
\begin{align}
u'(t) &= - u(t - 1/3) - u(t - 1/5), \qquad t \geq 0, \\
u(t) &= \delta_{0}(t), \qquad t \leq 0.
\end{align}
```

for ``t \in [0, 100]`` with the default fixed-point algorithm
and with Anderson acceleration.

The functions returns a named tuple of both solutions, called
`default` and `anderson`, respectively.

## References

Anderson, D. G. (1965) “Iterative Procedures for Nonlinear Integral Equations,” Journal of the ACM. Association for Computing Machinery (ACM). [doi: 10.1145/321296.321305](https://doi.org/10.1145/321296.321305).

Walker, H. F. and Ni, P. (2011) “Anderson Acceleration for Fixed-Point Iterations,” SIAM Journal on Numerical Analysis. Society for Industrial & Applied Mathematics (SIAM). [doi: 10.1137/10078356x](https://doi.org/10.1137/10078356x).
"""
function compare_anderson()
    # DDE problem
    prob = prob_dde_constant_2delays_long_scalar

    # Common algorithm and keyword arguments
    alg = Tsit5()
    kwargs = (reltol = 1e-3, abstol = 1e-6)

    # Solve the DDE with the default fixed-point iteration
    default_sol = solve(prob, MethodOfSteps(alg);
                        kwargs...)

    # Solve the DDE with Anderson acceleration
    anderson_sol = solve(prob,
                         MethodOfSteps(alg;
                                       fpsolve = NLAnderson());
                         kwargs...)

    return (; default = default_sol,
            anderson = anderson_sol)
end
