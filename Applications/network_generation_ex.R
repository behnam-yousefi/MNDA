library(mnda)
library(keras)

## Generate graph
myNet = network_gen(N_nodes = 100, N_var_nodes = 5, N_var_nei = 90, noise_sd = .01)
graphData = myNet[["data_graph"]]
edge_list = graphData[,1:2]
edge_weight = graphData[,3:4]

embeddingSpaceList = mnda_embedding_2layer(graphData, train.rep = 50, walk.rep=100)

Results = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "none")

myNet$var_nodes
Results$significant_nodes_index
Results$high_ranked_nodes_index
Results$high_var_nodes_index
