rm(list=ls())

library(mnda)

setwd("~/Desktop/R_Root/MNDA/usage_examples/")
graph_data = readRDS("Data/ISN_1.rds")
y = readRDS("Data/outcome_1.rds")

embeddingSpaceList = mnda_embedding(graph_data, train.rep=50, walk.rep=10,
                                           epochs = 20, batch.size = 10,
                                           random.walk=FALSE)

Results = mnda_node_detection_isn(embeddingSpaceList, stat.test = "t.test",
                                  p.adjust.method = "BH")

embeddingSpaceList = mnda_embedding(graph.data=graph_data, outcome=c(1,2), indv.index = c(1,1), train.rep=5, walk.rep=5)
Dist = mnda_node_distance(embeddingSpaceList)
