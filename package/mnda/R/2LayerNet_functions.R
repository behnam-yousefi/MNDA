## Pipeline for multiplex network differential analysis
## The aim of the current functions:
## Perform MNDA on two networks with different outcomes
## (i.e. case-control, time point1-time point2, etc.)
## and test whether the distances are significant.
## Example of the usage is provided in ???

# Null hypothesis (H0): The pairwise distances are not significantly larger than that of
## a random two-layer network.
## To have the null distribution, we generate two layers of permuted networks (from each layer)
## and embed the together as a 4-layer network.

#' Calculate the embedding space for a two layer multiplex network
#'
#' @param graph.data dataframe of the graph data containing edge list and edge weights.
#' column 1 and 2 consisting of the edge list (undirected).
#' column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
#' @param edge.threshold numeric value to set edge weights below the threshold to zero (default: 0). the greater edge weights do not change.
#' @param train.rep numeric value to set the number of EDNN training repeats (default: 50).
#' @param embedding.size the dimension of embedding space, equal to the number of the bottleneck hidden nodes (default: 5).
#' @param epochs maximum number of pocks. An early stopping callback with a patience of 5 has been set inside the function (default = 10).
#' @param batch.size batch size for learning (default = 5).
#' @param l2reg the coefficient of L2 regularization for the input layer (default = 0).
#'
#' @return a list of embedding spaces for each node.
#' @export
#'
#' @examples
#' embeddingSpaceList = mnda_embedding_2layer(graph.data)
#'
mnda_embedding_2layer = function(graph.data, edge.threshold=0, train.rep=50,
                       embedding.size=5, epochs=10, batch.size=5, l2reg=0){

  ### Step1) Read Graph Data ###
  ## The graph.data format is a data.frame:
  ## column 1 and 2 consisting of the edge list (undirected)
  ## column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
  EdgeList = data[,1:2]
  EdgeWeights = data[,3:4]
  outcome = colnames(EdgeWeights)

  ## generate the two null graphs by permuting edge weights of the original graphs:
  EdgeWeights = cbind(EdgeWeights,
                      sample(EdgeWeights[,1], nrow(EdgeWeights), replace = TRUE),
                      sample(EdgeWeights[,2], nrow(EdgeWeights), replace = TRUE))
  outcome = c(outcome, paste0(outcome, "_perm"))

  ### Step2) Preparing the input and output of the EDNN for all the graphs ###
  XY = ednn_IOprepare(edge.list = EdgeList, edge.weight = EdgeWeights,
                        outcome = outcome, edge.threshold = edge.threshold)
  X = XY["X"]
  Y = XY["Y"]

  ### Step3) EDNN training ###
  ## Process:
  ## 1. Train EDNN
  ## 2. Calculate the embedding space
  ## To have a stable measure, we repeat the process several times,
  ## which results in several distance vectors in the next step

  train.rep = 50
  embeddingSpaceList = list()
  for (rep in 1:train.rep)
    embeddingSpaceList[[rep]] = EDNN(X ,Y, Xtest = X,
                                     embedding.size, epochs, batch.size, l2reg)

  embeddingSpaceList[["outcome"]] = outcome_node
  return(embeddingSpaceList)
}

  ### Step4) Calculating the distance of node pairs in the embedding space ###
  ## To find the significantly varying nodes in the 2-layer-network, the distance between
  ## the corresponding nodes are calculated along with the null distribution.
  ## The null distribution is obtained based on the pairwise distances on null graphs.

#' Detecting the nodes whose local neighbors change bweteen the two conditions.
#'
#' @param embeddingSpaceList a list obtained by the \code{mnda_embedding_2layer()} function.
#'
#' @return the highly varibale nodes
#' @export
#'
#' @examples
#' embeddingSpaceList = mnda_embedding_2layer(graph.data)
#' Nodes = mnda_node_detection_2layer(embeddingSpaceList)
#'
mnda_node_detection_2layer = function(embeddingSpaceList){

  outcome_node = embeddingSpaceList[["outcome"]]
  outcome = unique(outcome_node)

  Rep = length(embeddingSpaceList)-1
  Dist = matrix(0, Rep, N_nodes)
  Dist_null = matrix(0, Rep, N_nodes)
  P_value = matrix(0, Rep, N_nodes)
  colnames(Dist) = names(V(graph))
  for (rep in 1:Rep){
    embeddingSpace = embeddingSpaceList[[rep]]
    embeddingSpace_1 = embeddingSpace[outcome_node==outcome[1], ]
    embeddingSpace_2 = embeddingSpace[outcome_node==outcome[2], ]
    for (i in 1:N_nodes)
      Dist[rep,i] = Distance(embeddingSpace_1[i,], embeddingSpace_2[i,], method = "cosine")

    ## Null distribution:
    embeddingSpace_1 = embeddingSpace[outcome_node==outcome[3], ]
    embeddingSpace_2 = embeddingSpace[outcome_node==outcome[4], ]
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
}
