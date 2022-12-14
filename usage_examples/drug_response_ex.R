rm(list=ls())

library(mnda)

setwd("~/Desktop/R_Root/MNDA/usage_examples/")
data = readRDS("Data/GCN2Layer_data_lung_tamoxifen_2000genes.rds")
X = data[[1]]
y = data[[2]]

adj_res = abs(cor(X[y=="res",]))
adj_nonres = abs(cor(X[y=="non_res",]))

diag(adj_res) = 0
diag(adj_nonres) = 0

adj_list = list(adj_res, adj_nonres)
graph_data = as.mnda.graph(adj_list, outcome = c("res","non_res"))

embeddingSpaceList = mnda_embedding_2layer(graph_data, edge.threshold = .1,
                                           train.rep = 50, walk.rep = 10,
                                           epochs = 20, batch.size = 10,
                                           random.walk = FALSE, null.perm = FALSE)

mnda_output = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "bonferroni", alpha = .05)

Nodes = mnda_output$high_var_nodes
# write.table(Nodes, "~/Desktop/nodes2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

subgraph_difference_plot(mnda.graph = graph_data, node.importance = mnda_output$rank_sum_dist,
                         n.var.nodes = 5, n.neigh = 10, diff.threshold = .1)
