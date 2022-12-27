rm(list=ls())

library(mnda)
setwd("~/Desktop/R_Root/MNDA/usage_examples/")

data = data.frame(readRDS("Data/ISN_net.rds"))
nodeList = t(sapply(rownames(data), function(x) strsplit(x,"_")[[1]]))

y = colnames(data)
y = t(data.frame(strsplit(y, "_")))

graph_data = cbind(nodeList, data)
embeddingSpaceList = mnda_embedding(graph_data, outcome = y, train.rep=5, walk.rep=10,
                                           epochs = 5, batch.size = 20,
                                           random.walk=FALSE)

Dist = mnda_node_distance(embeddingSpaceList)
