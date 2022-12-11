rm(list=ls())

library(mnda)

## Generate graph
myNet = network_gen(N_nodes = 100, N_var_nodes = 5, N_var_nei = 90, noise_sd = .01)
graph_data = myNet[["data_graph"]]

embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep = 50, walk.rep = 100,
                                           random.walk = TRUE, null.perm = TRUE)

mnda_output = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "none")

myNet$var_nodes
mnda_output$significant_nodes_index
mnda_output$high_ranked_nodes_index
mnda_output$high_var_nodes_index
hist(mnda_output$p_values)
