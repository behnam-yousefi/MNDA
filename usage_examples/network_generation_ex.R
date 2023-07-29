rm(list=ls())

library(mnda)

## Generate graph
myNet = network_gen(n.nodes = 100, n.var.nodes = 50, n.var.nei = 90, noise.sd = .01)
graph_data = myNet[["data_graph"]]

embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep = 10, walk.rep = 100,
                                           random.walk = TRUE, null.perm = TRUE, verbose = TRUE)

mnda_output = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "none")

myNet$var_nodes
mnda_output$significant_nodes_index
mnda_output$high_ranked_nodes_index
mnda_output$high_var_nodes_index
hist(mnda_output$p_values)
