rm(list=ls())

library(mnda)
setwd("~/Desktop/R_Root/MNDA/usage_examples/")

data = readRDS("Data/ISN_BCG.rds")
nodeList = t(sapply(rownames(data), function(x) strsplit(x,"_")[[1]]))

y = colnames(data)
y = data.frame(t(data.frame(strsplit(y, "_"))))
colnames(y) = c("ID", "Stim", "Sex")
sex = y[duplicated(y$ID), "Sex"]

graph_data = cbind(nodeList, data)
embeddingSpaceList = mnda_embedding(graph_data, outcome = y$Stim, indv.index = y$ID,
                                    train.rep=2, walk.rep=10, epochs=5, batch.size=100,
                                    random.walk=FALSE)

Dist = mnda_node_distance(embeddingSpaceList)
Pval = mnda_distance_test_isn(Dist, sex)



