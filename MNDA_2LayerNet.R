## Multiplex Graph Representation Learning and
## multiplex network differential analysis
## By: Behnam Yousefi

## The aim of the current script:
## perform MNDA on two networks with different outcomes 
## (i.e. case-control, time point1-time point2, etc.)

rm(list = ls())
setwd("~/Desktop/R_Root/MNDA/")
source('MNDA_functions.R')
source('Figure_subgraph_func.R')

library(igraph)
library(MASS)

### Step1) Read Graph Data ###
# The graph format is a data.frame:
# column 1 and 2 consisting of the edge list (undirected)
# column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
outcome = c(1,0)                # 1: responder, 0: nonresponder
data = readRDS("Data/MNDA-drug/CD_TNF_w14_Global.rds")
NodeList = data[,1:2]
EdgeWeights = data[,3:4]
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

Rep = 10
embeddingSpaceList = list()
for (rep in 1:Rep)
  embeddingSpaceList[[rep]] = EDNN(X ,Y, Xtest = X, latentSize = 5, 
                                   epochs = 20, batch_size = 5, l2reg = .0001)

# Plot 
embeddingSpace = embeddingSpaceList[[1]]
plot(embeddingSpace[,1:2], pch = 20, cex = .5, col = outcome_node+1)

# embeddingSpaceList[["outcome"]] = outcome_node
# saveRDS(embeddingSpaceList, "Data/Embedding_Space/Embedding_Space_1.rds")
# embeddingSpaceList = readRDS("Data/Embedding_Space/Embedding_Space_1.rds")

### Step4) Calculating the distance of node pairs in the embedding space ###
## To find the highly variable nodes in a 2-layer-network, the distance between the corresponding
## nodes between the two layers is calculated.
Dist = matrix(0, Rep, N_nodes)
colnames(Dist) = names(V(graph))
for (rep in 1:Rep){
  embeddingSpace = embeddingSpaceList[[rep]]
  embeddingSpace_1 = embeddingSpace[outcome_node==1, ]
  embeddingSpace_2 = embeddingSpace[outcome_node==0, ]
  for (i in 1:N_nodes)
    Dist[rep,i] = Distance(embeddingSpace_1[i,], embeddingSpace_2[i,], method = "cosine")
}

hist(Dist[,], 20)

### Step5) Find the highly and lowly variant nodes using
### a rank sum-based method

Rank_dist = matrix(0, Rep, N_nodes)
colnames(Rank_dist) = names(V(graph))
for (rep in 1:Rep)
  Rank_dist[rep,] = Rank(Dist[rep,], decreasing = FALSE)

Rank_sum_dist = apply(Rank_dist, 2, sum)

plot(sort(Rank_sum_dist), pch = 20)
abline(h = 3810, col = "red")
high_var_nodes = order(Rank_sum_dist, decreasing = TRUE)[1:9]

print(high_var_nodes)
print(names(Rank_sum_dist)[high_var_nodes])
Node_set = names(Rank_sum_dist)[high_var_nodes]

a = c(131, 136, 130, 40, 134, 124, 99, 112, 104)
b = c(131, 130, 40, 124, 115, 136, 104, 134, 135)

### Step6) Plot subgraph limited to the high-var nodes ###
W = as.numeric(EdgeWeights[,1]) - as.numeric(EdgeWeights[,2]) 
Threshold = 0
graph_plot = simplify(set.edge.attribute(graph, "weight", index=E(graph), W))
graph_plot = delete_edges(graph_plot, E(graph_plot)[abs(E(graph_plot)$weight) < Threshold])
Labels = names(Rank_sum_dist)
names(Labels) = names(Rank_sum_dist)

subgraph_plot(graph_plot, Node_set, Labels)

