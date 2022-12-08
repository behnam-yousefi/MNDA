rm(list=ls())

library(mnda)

## Generate graph
myNet = network_gen(N_nodes = 100, N_var_nodes = 5, N_var_nei = 90, noise_sd = .01)
graph_data = myNet[["data_graph"]]

embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep=50, walk.rep=100,
                                           random.walk=TRUE, null.perm = TRUE)

Results = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "none")

myNet$var_nodes
Results$significant_nodes_index
Results$high_ranked_nodes_index
Results$high_var_nodes_index
