alg_order(alg::AbstractMethodOfStepsAlgorithm) = alg_order(alg.alg)
isconstrained{constrained}(alg::AbstractMethodOfStepsAlgorithm{constrained}) = constrained
