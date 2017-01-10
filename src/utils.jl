#=
lags = [1//5,1//3,1//2]
n = length(lags)
order = 3
start_val = 1.05
function compute_discontinuity_tree(lags,order,start_val)
  n = length(lags)
  locs = vcat([vcat(unique.(permutations.([i for i in with_replacement_combinations((0:n),n) if sum(i) == j]))...) for j in 1:order]...)
  vals = unique([start_val + sum(loc.*lags) for loc in locs])
end
compute_discontinuity_tree(lags,order,start_val)
=#

function compute_discontinuity_tree(lags,alg,start_val)
  n = length(lags)
  locs = vcat([vcat(unique.(permutations.([i for i in with_replacement_combinations((0:alg_order(alg)),n) if sum(i) == j]))...) for j in 1:alg_order(alg)]...)
  vals = unique([start_val + sum(loc.*lags) for loc in locs])
end
