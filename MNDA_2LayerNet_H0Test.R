## Multiplex Graph Representation Learning and
## multiplex network differential analysis
## By: Behnam Yousefi

## The aim of the current script:
## Perform MNDA on two networks with different outcomes 
## (i.e. case-control, time point1-time point2, etc.)
## and test whether the distances are significant.
#
# Null hypothesis (H0): The pairwise distances are not significantly larger than that of
## a random two-layer network.
## To have the null distribution, we generate two layers of permuted networks (from each layer)
## and embed the together as a 4-layer network.

rm(list = ls())
setwd("~/Desktop/R_Root/MNDA/")
source('MNDA_functions.R')
source('Figure_subgraph_func.R')

library(igraph)
library(MASS)
library(aggregation)

### Step1) Read Graph Data ###
## The graph format is a data.frame:
## column 1 and 2 consisting of the edge list (undirected)
## column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
outcome = c(1,0,2,3) # 1: responder, 0: nonresponder, 2: permuted resp., 3: perm. nonr.
data = readRDS("Data/MNDA-drug/CD_TNF_w14_Global.rds")
NodeList = data[,1:2]
EdgeWeights = data[,3:4]

## generate the two null graphs by permuting edge weights of the original graphs:
EdgeWeights = cbind(EdgeWeights, 
                    sample(EdgeWeights[,1], nrow(EdgeWeights), replace = TRUE),
                    sample(EdgeWeights[,2], nrow(EdgeWeights), replace = TRUE))

graph = graph(t(NodeList), directed = FALSE)

N_nodes = length(V(graph))
N_graph = ncol(EdgeWeights)

### Step2) Preparing the input and output of the EDNN for all the graphs ###
Threshold = 0
X = c()
Y = c()
outcome_node = c()

for (i in 1:N_graph){
  ## Set layer-specific weights for each graph and
  ## [optional] perform a thresholding if needed
  W_i = as.numeric(EdgeWeights[,i])
  graph_i = simplify(set.edge.attribute(graph, "weight", index=E(graph), W_i))
  graph_i = delete_edges(graph_i, E(graph_i)[E(graph_i)$weight < Threshold])
  
  ## Step2.1) Input: Adjacency matrix calculation
  Adj_i = as.matrix(as_adj(graph_i,  attr = "weight"))

  ## Step2.2) Output: Perform the fixed-length random walk 
  ## and calculating the node visit probabilities.
  # Two options exist:
  # 1.repetitive simple random walks
  # 2.repetitive weighted random walks (specific to weighted graphs)
  RW = RepRandomWalk (graph_i, Nrep = 100, Nstep = 5, weighted_walk = TRUE)

  ## Step2.3) Make it multilayer for EDNN
  X = rbind(X, Adj_i)
  Y = rbind(Y, RW$Probabilities)
  
  outcome_node = c(outcome_node, rep(outcome[i], N_nodes))
}
X = X / (apply(X, 1, sum) + .000000001)

### Step3) EDNN training ###
## Process:
## 1. Train EDNN
## 2. Calculate the embedding space
## To have a stable measure, we repeat the process several times,
## which results in several distance vectors in the next step

Rep = 50
embeddingSpaceList = list()
for (rep in 1:Rep)
  embeddingSpaceList[[rep]] = EDNN(X ,Y, Xtest = X, latentSize = 5, 
                                   epochs = 10, batch_size = 5, l2reg = .0001)

# Plot 
embeddingSpace = embeddingSpaceList[[1]]
plot(embeddingSpace[,4:5], pch = 20, cex = .5, col = outcome_node+1)

# embeddingSpaceList[["outcome"]] = outcome_node
# saveRDS(embeddingSpaceList, "Data/Embedding_Space/Embedding_Space_2.rds")
# embeddingSpaceList = readRDS("Data/Embedding_Space/Embedding_Space_2.rds")

### Step4) Calculating the distance of node pairs in the embedding space ###
## To find the significantly varying nodes in the 2-layer-network, the distance between 
## the corresponding nodes are calculated along with the null distribution.
## The null distribution is obtained based on the pairwise distances on null graphs.
Dist = matrix(0, Rep, N_nodes)
Dist_null = matrix(0, Rep, N_nodes)
P_value = matrix(0, Rep, N_nodes)
colnames(Dist) = names(V(graph))
for (rep in 1:Rep){
  embeddingSpace = embeddingSpaceList[[rep]]
  embeddingSpace_1 = embeddingSpace[outcome_node==1, ]
  embeddingSpace_2 = embeddingSpace[outcome_node==0, ]
  for (i in 1:N_nodes)
    Dist[rep,i] = Distance(embeddingSpace_1[i,], embeddingSpace_2[i,], method = "cosine")
  
  ## Null distribution:
  embeddingSpace_1 = embeddingSpace[outcome_node==2, ]
  embeddingSpace_2 = embeddingSpace[outcome_node==3, ]
  for (i in 1:N_nodes)
    Dist_null[rep,i] = Distance(embeddingSpace_1[i,], embeddingSpace_2[i,], method = "cosine")
  
  for (i in 1:N_nodes)
    P_value[rep,i] = wilcox.test(Dist_null[rep,], y = Dist[rep,i], 
                                 alternative = "less")$p.value
}

hist(Dist[,], 20)
hist(P_value[,6], 50, xlim = c(.01,1))
abline(v = .05, col = "red")

P_value_aggr = apply(P_value, 2, fisher)
hist(P_value_aggr, 50)
abline(v = .05, col = "red")
which(P_value_aggr<.05)


### Step5) Find the highly and lowly variant nodes using
### a rank sum-based method
Rank_dist = matrix(0, Rep, N_nodes)
colnames(Rank_dist) = names(V(graph))
for (rep in 1:Rep)
  Rank_dist[rep,] = Rank(Dist[rep,], decreasing = FALSE)

Rank_sum_dist = apply(Rank_dist, 2, sum)

plot(Rank_sum_dist, -log10(P_value_aggr))
abline(h = -log10(.01), col = "red")
abline(v = 6000, col = "red")

plot(sort(Rank_sum_dist), pch = 20)
abline(h = 3810, col = "red")
high_var_nodes = order(Rank_sum_dist, decreasing = TRUE)[1:5]
low_var_nodes = order(Rank_sum_dist, decreasing = FALSE)[1:11]

print(high_var_nodes)
print(names(V(graph))[high_var_nodes])
Node_set = names(V(graph))[high_var_nodes]

### Step6) Plot subgraph limited to the high-var nodes ###
W = as.numeric(EdgeWeights[,1]) - as.numeric(EdgeWeights[,2]) 
Threshold = 0
graph_plot = simplify(set.edge.attribute(graph, "weight", index=E(graph), W))
graph_plot = delete_edges(graph_plot, E(graph_plot)[abs(E(graph_plot)$weight) < Threshold])
Labels = names(V(graph))

subgraph_plot(graph_plot, Node_set, Labels, font_size = 4, margin = 2.5)

