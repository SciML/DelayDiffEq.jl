alg_order(alg::AbstractMethodOfStepsAlgorithm) = alg_order(alg.alg)

isconstrained(alg::AbstractMethodOfStepsAlgorithm{constrained}) where constrained =
    constrained
