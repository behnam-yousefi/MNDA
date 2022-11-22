rm(list=ls())

library(mnda)
library(keras)
library(aggregation)

setwd("~/Desktop/R_Root/MNDA/Applications/")
data = readRDS("Data/graph_data_lung_tamoxifen.rds")
X = data[[1]]
y = data[[2]]

adj_res = abs(cor(X[y=="res",]))
adj_nonres = abs(cor(X[y=="non_res",]))

adj_list = list(adj_res, adj_nonres)
graph_data = as.mnda.graph(adj_list, outcome = c("res","non_res"))

embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep=10, walk.rep=10,
                                           random.walk=TRUE, null.perm = TRUE)

Results = mnda_node_detection_2layer(emb_list, p.adjust.method = "none")

