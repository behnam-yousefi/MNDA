rm(list=ls())

library(mnda)

setwd("~/Desktop/R_Root/MNDA/usage_examples/")
graph_data = readRDS("Data/ISN_1.rds")
y = readRDS("Data/outcome_1.rds")

embeddingSpaceList = mnda_embedding(graph_data, outcome = y,train.rep=5, walk.rep=10,
                                           epochs = 5, batch.size = 20,
                                           random.walk=FALSE)

Dist = mnda_node_distance(embeddingSpaceList)
