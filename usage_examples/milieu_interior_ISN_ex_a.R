rm(list=ls())

library(mnda)
setwd("~/Desktop/R_Root/MNDA/usage_examples/")

## Read ISN data and construct the node list
data = readRDS("Data/ISN_BCG.rds")
nodeList = t(sapply(rownames(data), function(x) strsplit(x,"_")[[1]]))

## Construct the phenotype data
y = colnames(data)
y = data.frame(t(data.frame(strsplit(y, "_"))))
colnames(y) = c("ID", "Stim", "Sex")

## Construct the two aggregated networks
data_agg = cbind(apply(data[,y$Stim=="Null"], 1, mean),
                 apply(data[,y$Stim=="BCG"], 1, mean))

## Run the algorithm
## 1. Embed nodes
graph_data = cbind(nodeList, data_agg)
embeddingSpaceList = mnda_embedding_2layer(graph_data, edge.threshold = .1,
                                           train.rep = 50, epochs = 25, batch.size = 10,
                                           random.walk = FALSE, null.perm = FALSE)

# saveRDS(embeddingSpaceList, file = "Data/Embeddings/embeddingSpaceList_BCG_a.rds")
# embeddingSpaceList = readRDS("Data/Embeddings/embeddingSpaceList_S.Aureus.rds")

## 2. Calculate distances and p.values
mnda_output = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "bonferroni", alpha = .01)

# Nodes = mnda_output$high_var_nodes
# write.table(Nodes, "Data/nodes_drug.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Plot the difference sub-graph
plt = subgraph_difference_plot(mnda.graph = graph_data, node.importance = mnda_output$rank_sum_dist,
                               n.var.nodes = 10, n.neigh = 10, diff.threshold = .2, edge.width = 3)
plt

# pdf("../Figures/subgraph_BCG_a.pdf")
# plt
# dev.off()
