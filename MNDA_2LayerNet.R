## Multiplex Graph Representation Learning and
## multiplex network differential analysis
## By: Behnam Yousefi

## The aim of the current script:
## perform MNDA on two networks with different outcomes 
## (i.e. case-control, time point1-time point2, etc.)

rm(list = ls())
setwd("~/Desktop/R_Root/MNDA/")
source('~/Desktop/R_Root/MNDA/MNDA_functions.R')

library(igraph)
library(MASS)
library(data.table)

### Step1) Read Graph Data ###
# The graph format is a data.frame:
# column 1 and 2 consisting of the edge list (undirected)
# column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
outcome = a
data = data.frame(fread(file="Data/Resulting_net_notNULL_MAGMACONF6M.txt",sep = " "))
NodeList = data[,1:3]
EdgeWeights = data[,3:4]
graph = graph(t(NodeList), directed = FALSE)

N_nodes = length(V(graph1))
N_graph = ncol(EdgeWeights)

### Step2) Preparing the input and output of the EDNN for all the graphs ###
Threshold = 0
X = c()
Y = c()
outcome_node = c()

for (i in N_graph){
  ## Set layer-specific weights for each graph and
  ## [optional] perform a thresholding if needed
  W_i = EdgeWeights[,i]
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
  Y = rbind(Y, RW$P)
  
  outcome_node = c(outcome_node, rep(outcome[i], nrow(N_nodes)))
}
X = X / (apply(X, 1, sum) + .000000001)

colnames(X) = paste0("V",as.character(1:ncol(X)))
colnames(Y) = paste0("V",as.character(1:ncol(Y)))

### Step3) EDNN training and calculating the distance of node pairs ###
## Process:
## 1. Train EDNN
## 2. Calculate the embedding space
## 3. Calculate the distance between the node pairs
## To have a stable measure, we repeat the process several times,
## which results in several distance vectors

Rep = 30
Dist = matrix(0, Rep, N)
for (rep in 1:Rep){
  # latentSpace = EDNN(X, Y, Xtest = X, latentSize = 2, epochs = 10, batch_size = 100)
  latentSpace = EDNN(X ,Y, Xtest = X, latentSize = 10, epochs = 50, batch_size = 5, l2reg = .0001)
  
  # plot(latentSpace[,1:2],pch = 20, cex = .1)
  
  latentSpace_1 = latentSpace[outcome_node==1, ]
  latentSpace_2 = latentSpace[outcome_node==2, ]
  for (i in 1:N)
    Dist[rep,i] = Distance (latentSpace_1[i,], latentSpace_2[i,], method = "cosine")
}

### Step4) Find the highly and lowly variant nodes using
### a rank sum-based method

Rank_Dist = matrix(0, Rep, N)
for (rep in 1:Rep){
  Dist[rep,]
}


