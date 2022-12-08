rm(list=ls())

library(mnda)

setwd("~/Desktop/R_Root/MNDA/usage_examples/")
data = readRDS("Data/GCN2Layer_data_lung_tamoxifen_1000genes.rds")
X = data[[1]]
y = data[[2]]

adj_res = abs(cor(X[y=="res",]))
adj_nonres = abs(cor(X[y=="non_res",]))

adj_list = list(adj_res, adj_nonres)
graph_data = as.mnda.graph(adj_list, outcome = c("res","non_res"))

embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep=50, walk.rep=10,
                                           epochs = 20, batch.size = 10,
                                           random.walk=FALSE, null.perm = FALSE)

Results = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "BH")
Nodes = Results$high_var_nodes
# write.table(Nodes, "~/Desktop/nodes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

NodeSet = c("IL1RN", "CCND1", "HES1", "PTMS", "ACOT4", "KLF5")

graph_to_plot = cbind(graph_data[,1:2], W = graph_data[,3] - graph_data[,4])
hist(graph_to_plot$W)

G = mnda::as.igraph(graph_to_plot, 1.15)

pdf("~/Desktop/cancer.pdf")
subgraph_plot(G, NodeSet)
dev.off()
