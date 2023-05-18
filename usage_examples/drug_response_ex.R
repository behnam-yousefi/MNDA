rm(list=ls())

library(mnda)
setwd("~/Desktop/R_Root/MNDA/usage_examples/")

## Read graph data
data = readRDS("Data/GCN2Layer_data_lung_tamoxifen_2000genes.rds")
X = data[[1]]
y = data[[2]]

## Construct the phenotype data
adjRes = abs(cor(X[y=="res",]))
adjNonres = abs(cor(X[y=="non_res",]))

## Construct adjacency matrices of gene coexpression network
diag(adjRes) = 0
diag(adjNonres) = 0

## Construct multiplex network
adjList = list(adjRes, adjNonres)
graphData = as_mnda_graph(adjList, outcome = c("res","non_res"))

## Run the algorithm
## 1. Embed nodes
embeddingSpaceList = mnda_embedding_2layer(graphData, edge.threshold = .1,
                                           train.rep = 50, epochs = 25, batch.size = 10,
                                           random.walk = FALSE, null.perm = FALSE)

# saveRDS(embeddingSpaceList, file = "Data/Embeddings/embeddingSpaceList_drug.rds")
embeddingSpaceList = readRDS("Data/Embeddings/embeddingSpaceList_drug.rds")

## 2. Calculate distances and p.values
mnda_output = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "bonferroni", alpha = .01)

# Nodes = mnda_output$high_var_nodes
# write.table(Nodes, "Data/nodes_drug.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Plot the difference sub-graph
plt = subgraph_difference_plot(mnda.graph = graphData, node.importance = mnda_output$rank_sum_dist,
                         n.var.nodes = 10, n.neigh = 10, diff.threshold = .2, edge.width = 3)
plt

# pdf("../Figures/subgraph_drug.pdf")
# plt
# dev.off()

## Expression test
pValExpr = apply(X, 2, FUN = function(x){
  p.val = t.test(x[y=="res"], x[y=="non_res"])$p.value
  return(p.val)
})

hist(pValExpr,100)
