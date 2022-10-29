## Multiplex Graph Representation Learning and
## multiplex netwoek differential analysis
## Graph Representation Learning
## By: Behnam Yousefi

rm(list = ls())
setwd("~/Desktop/R_Root/MNDA/")
source('~/Desktop/R_Root/MNDA/MNDA_functions.R')

library(igraph)
library(MASS)
library(data.table)

### Step1) Read Graph Data ###
# Two graph data is needed for each set (case-control, time point1-time point2, etc.)
# The graph format is a node list as a two-row matrix;
# [optional] the edge weights can be added as the 3rd row to the imput matrix.
# Graph 1
data = data.frame(fread(file="Data/Resulting_net_notNULL_MAGMACONF6M.txt",sep = " " ))
NodeList = graph[1:2,]
EdgeWeight = graph[3,]
graph1 = graph(NodeList, directed = FALSE)
graph1 = simplify(set.edge.attribute(graph1, "weight", index=E(graph1), EdgeWeight))

# Graph 2
data = data.frame(fread(file="Data/Resulting_net_notNULL_MAGMACONF6M.txt",sep = " " ))
NodeList = graph[1:2,]
EdgeWeight = graph[3,]
graph2 = graph(NodeList, directed = FALSE)
graph2 = simplify(set.edge.attribute(graph2, "weight", index=E(graph2), EdgeWeight)) 

assertthat::assert_that(length(V(graph1)) == length(V(graph2)))
N_nodes = length(V(graph1))

# Graph weight thresholding if needed
Threshold = 0
graph1 = delete_edges(graph1, E(graph1)[E(graph1)$weight < Threshold])
graph2 = delete_edges(graph2, E(graph2)[E(graph2)$weight < Threshold])

### Step2) Preparing the input and output of the EDNN  ###
## Step2.1) Input: Adjacency matrix calculation
A1 = as.matrix(as_adj(graph1,  attr = "weight"))
A2 = as.matrix(as_adj(graph2,  attr = "weight"))

## Step2.2) Output: Perform the fixed-length random walk 
## and calculating the node visit probabilities.
# Two options exist:
# 1.repetitive simple random walks
# 2.repetitive weighted random walks (specific to weighted graphs)
RW1 = RepRandomWalk (graph1, Nrep = 100, Nstep = 5, weighted_walk = TRUE)
RW2 = RepRandomWalk (graph2, Nrep = 100, Nstep = 5, weighted_walk = TRUE)

# Summary
X = rbind(A1, A2)
X = X / (apply(X, 1, sum) + .000000001)
Y = rbind(RW1$P, RW2$P)

colnames(X) = paste0("V",as.character(1:ncol(X)))
colnames(Y) = paste0("V",as.character(1:ncol(Y)))

individual_no = rep(1:nrow(X)/2, 2)
set = c(rep(1, nrow(X)/2), rep(2, nrow(X)/2))

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
  
  latentSpace_1 = latentSpace[set==1, ]
  latentSpace_2 = latentSpace[set==2, ]
  for (i in 1:N)
    Dist[rep,i] = Distance (latentSpace_1[i,], latentSpace_2[i,], method = "cosine")
}

### Step4) Find the highly and lowly variant nodes using
### a rank sum-based method

Rank_Dist = matrix(0, Rep, N)
for (rep in 1:Rep){
  Dist[rep,]
}


