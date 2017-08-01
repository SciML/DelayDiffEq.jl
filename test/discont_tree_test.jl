using DelayDiffEq, OrdinaryDiffEq, DataStructures, Base.Test

lags = [1//5, 1//2]
alg = BS3()
start_val = 1

disc_tree = sort(DelayDiffEq.compute_discontinuity_tree(lags, alg, start_val))

true_val = sort(1 .+ [1//2, 1//5, 1//1, 2//5, 7//10, 3//2, 3//5, 6//5, 9//10])
@test disc_tree == true_val
