rm(list=ls())

library(mnda)
setwd("~/Desktop/R_Root/MNDA/usage_examples/")

data = readRDS("Data/ISN_S.Aureus.rds")
nodeList = t(sapply(rownames(data), function(x) strsplit(x,"_")[[1]]))

y = colnames(data)
y = data.frame(t(data.frame(strsplit(y, "_"))))
colnames(y) = c("ID", "Stim", "Sex")
sex = y[duplicated(y$ID), "Sex"]

graph_data = cbind(nodeList, data)
embeddingSpaceList = mnda_embedding(graph_data, outcome = y$Stim, indv.index = y$ID,
                                    train.rep=50, walk.rep=10, epochs=10, batch.size=50,
                                    random.walk=FALSE)

# saveRDS(embeddingSpaceList, file = "Data/Embeddings/embeddingSpaceList_S.Aureus.rds")
# embeddingSpaceList = readRDS("Data/Embeddings/embeddingSpaceList_S.Aureus.rds")

Dist = mnda_node_distance(embeddingSpaceList)
Pval = mnda_distance_test_isn(Dist, sex, p.adjust.method = "bonferroni")
# sum(Pval<.01)
# Pval = sort(Pval, decreasing = FALSE)
# TopVarGenes = names(Pval[1:10])

graph_to_plot = cbind(nodeList,
                      apply(data[,y$Stim == "Null"], 1, mean),
                      apply(data[,y$Stim != "Null"], 1, mean))

subgraph_difference_plot(mnda.graph = graph_to_plot, node.importance = -log(Pval),
                         n.var.nodes = 10, n.neigh = 10, diff.threshold = .1, edge.width = 5)



