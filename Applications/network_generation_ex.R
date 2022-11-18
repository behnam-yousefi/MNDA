library(mnda)
library(keras)

## Generate graph
myNet = network_gen(N_nodes = 100, N_var_nodes = 5, N_var_nei = 90, noise_sd = .01)
graphData = myNet[["data_graph"]]
edge_list = graphData[,1:2]
edge_weight = graphData[,3:4]

embeddingList = mnda_embedding_2layer(graphData, train.rep = 2, walk.rep=10)
