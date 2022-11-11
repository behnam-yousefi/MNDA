## Multiplex Network Simulation
## In this script we generate random pairs of gene co-expression networks, which
## are different only in a few (pre-set) number of nodes.

rm(list = ls())
setwd("~/Desktop/R_Root/MNDA/")

N_nodes = 100
N_var_nodes = 5
noise_sd = .1
N_var_nei = 90

Adj1 = matrix(runif(N_nodes^2), N_nodes)
Adj2 = Adj1

var_node_set = sample(N_nodes, N_var_nodes, replace = FALSE)
for (node in var_node_set){
  var_nei_set = sample((1:N_nodes)[-node], N_var_nei, replace = FALSE)
  new_edge_weights = runif(N_var_nei)
  Adj2[node, var_nei_set] = new_edge_weights
  Adj2[var_nei_set, node] = Adj2[node, var_nei_set]
}

noise = abs(matrix(rnorm(N_nodes^2, 0, noise_sd), N_nodes))
Adj2 = Adj2 + noise

diag(Adj1) = 0
diag(Adj2) = 0
Adj1[lower.tri(Adj1)] = t(Adj1)[lower.tri(Adj1)]
Adj2[lower.tri(Adj2)] = t(Adj2)[lower.tri(Adj2)]


NodeList = data.frame(V1 = rep(1:N_nodes, times = N_nodes),
                      V2 = rep(1:N_nodes, each = N_nodes))
EdgeWeights = data.frame(W1 = as.numeric(Adj1),
                         W2 = as.numeric(Adj2))
print(var_node_set)
