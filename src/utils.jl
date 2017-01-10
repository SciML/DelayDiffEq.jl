function compute_discontinuity_tree(lags, alg, start_val)
           start_val + unique(vcat((sum.(collect(with_replacement_combinations(lags, i))) for i in 1:alg_order(alg))...))
       end
