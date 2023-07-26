rm(list=ls())

library(mnda)
setwd("~/Desktop/R_Root/MNDA/usage_examples/")

## Read ISN data and construct the node list
data = readRDS("Data/ISN_S.Aureus.rds")
nodeList = t(sapply(rownames(data), function(x) strsplit(x,"_")[[1]]))

## Construct the phenotype data
y = colnames(data)
y = data.frame(t(data.frame(strsplit(y, "_"))))
colnames(y) = c("ID", "Stim", "Sex")
sex = y[duplicated(y$ID), "Sex"]

## Run the algorithm
## 1. Embed nodes
graph_data = cbind(nodeList, data)
embeddingSpaceList = mnda_embedding(graph_data, outcome = y$Stim, indv.index = y$ID,
                                    train.rep=5, walk.rep=1, epochs=1, batch.size=50,
                                    random.walk=FALSE)

# saveRDS(embeddingSpaceList, file = "Data/Embeddings/embeddingSpaceList_S.Aureus.rds")
embeddingSpaceList = readRDS("Data/Embeddings/embeddingSpaceList_S.Aureus.rds")

## 2. Calculate distances
Dist = mnda_node_distance(embeddingSpaceList)

## 3. Calculate p.valus
Pval = mnda_distance_test_isn(Dist, sex, p.adjust.method = "bonferroni")
# sum(Pval<.01)
# Pval = sort(Pval, decreasing = FALSE)
# TopVarGenes = names(Pval[1:10])

## Plot the difference sub-graph (difference between the average network of each condition)
graph_to_plot = cbind(nodeList,
                      apply(data[,y$Stim == "Null"], 1, mean),
                      apply(data[,y$Stim != "Null"], 1, mean))


plt = subgraph_difference_plot(mnda.graph = graph_to_plot, node.importance = -log(Pval),
                               n.var.nodes=10, n.neigh=10, edge.width=c(.5,1))
plt

# pdf("../Figures/subgraph_S.Aureus.pdf")
# plt
# dev.off()
