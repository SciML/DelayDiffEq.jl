alg_maximum_order(alg::AbstractMethodOfStepsAlgorithm) = alg_maximum_order(alg.alg)

isconstrained(alg::AbstractMethodOfStepsAlgorithm{constrained}) where constrained =
    constrained
