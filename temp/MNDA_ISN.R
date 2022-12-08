## Multiplex Graph Representation Learning and
## multiplex network differential analysis
## By: Behnam Yousefi

## The aim of the current script:
## perform MNDA on a set of ISNs (e.g. LIONESS - Kuijjer et al. 2019)
## Two conditions we have:
## 1) Paired ISNs: for each condition we have the same individuals (e.g. time-based analysis)
## 2) Unpaired ISNs: individuals are different between conditions (e.g. case-control)
##  To find the variable nodes in the Paired ISNs case, we can calculate distances between 
## corresponding nodes and compare which results in a set of distances for each node.
## A Wilcoxon test can be then be used to find the significant nodes.
##  For the Unpaired case, such a distance cannot be calculated. The only method in this case
## would be to use the embedding space of the nodes as input features of a predictive ML model
## y = wX, X [individual by features], |features| = N_nodes by N_dim. A leave-one-out CV can
## be used to train and test the model. Since this analysis is context dependent, we do not 
## provide any standard pipeline for the Unpaired ISNs case.

rm(list = ls())
setwd("~/Desktop/R_Root/MNDA/")
source('MNDA_functions.R')
source('Figure_subgraph_func.R')

library(igraph)
library(MASS)

### Step1) Read Graph Data ###
## The graph format is a data.frame:
## column 1 and 2 consisting of the edge list (undirected)
## column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
outcome = readRDS("Data/MNDA-drug/CD_TNF_w14_outcome.rds")$response
data = readRDS("Data/MNDA-drug/CD_TNF_w14_ISN.rds")
NodeList = data[,1:2]
EdgeWeights = data[,3:ncol(data)]
graph = graph(t(NodeList), directed = FALSE)

N_nodes = length(V(graph))
N_graph = ncol(EdgeWeights)

### Step2) Preparing the input and output of the EDNN for all the graphs ###
Threshold = 0
X = c()
Y = c()
outcome_node = c()
individual_node = c()

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
  individual_node = c(individual_node, rep(i, N_nodes))
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
                                   epochs = 20, batch_size = 5, l2reg = 0)

# Plot 
embeddingSpace = embeddingSpaceList[[1]]
plot(embeddingSpace[,1:2], pch = 20, cex = .5, col = outcome_node+1)

# embeddingSpaceList[["outcome"]] = outcome_node
# saveRDS(embeddingSpaceList, "Data/Embedding_Space/Embedding_Space_1.rds")
# embeddingSpaceList = readRDS("Data/Embedding_Space/Embedding_Space_1.rds")
