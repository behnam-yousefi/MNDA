## Pipeline for multiplex network differential analysis
## The aim of the current functions:
## Perform MNDA on two networks with different outcomes
## (e.g. case-control, time point1-time point2, etc.)
## and test whether the distances are significant.
## Example of the usage is provided in the GitHub

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
#' @param walk.rep number of repeats for the random walk (default: 100).
#' @param n.steps number of the random walk steps (default: 5).
#' @param random.walk boolean value to enable the random walk algorithm (default: TRUE).
#' @param null.perm boolean to enable permuting two random graphs and embed them, along with the main two graphs, for the null distribution (default: TRUE).
#' @param demo a boolean vector to indicate this is a demo example or not
#' @param verbose if \emph{TRUE} a progress bar is shown.
#'
#' @return a list of embedding spaces for each node.
#' @export
#'
#' @examples
#' myNet = network_gen(N_nodes = 50, N_var_nodes = 5, N_var_nei = 40, noise_sd = .01)
#' graph_data = myNet[["data_graph"]]
#' embeddingSpaceList = mnda_embedding_2layer(graph.data=graph_data, train.rep=5, walk.rep=5)
#'
mnda_embedding_2layer = function(graph.data, edge.threshold=0, train.rep=50,
                       embedding.size=2, epochs=10, batch.size=5, l2reg=0,
                       walk.rep=100, n.steps=5, random.walk=TRUE, null.perm=TRUE, demo = TRUE, verbose=TRUE){

  assertthat::assert_that(train.rep >= 2)

  ### Step1) Read Graph Data ###
  ## The graph.data format is a data.frame:
  ## column 1 and 2 consisting of the edge list (undirected)
  ## column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
  EdgeList = graph.data[,1:2]
  EdgeWeights = graph.data[,3:4]
  outcome = colnames(EdgeWeights)

  ## generate the two null graphs by permuting edge weights of the original graphs:
  if (null.perm){
    EdgeWeights = cbind(EdgeWeights,
                        sample(EdgeWeights[,1], nrow(EdgeWeights), replace = TRUE),
                        sample(EdgeWeights[,2], nrow(EdgeWeights), replace = TRUE))
    outcome = c(outcome, paste0(outcome, "_perm"))
  }

  ### Step2) Preparing the input and output of the EDNN for all the graphs ###
  XY = ednn_IOprepare(edge.list = EdgeList, edge.weight = EdgeWeights,
                      walk.rep = walk.rep, n.steps = n.steps, random.walk = random.walk,
                      outcome = outcome, edge.threshold = edge.threshold, verbose = verbose)
  X = XY[["X"]]
  Y = XY[["Y"]]
  outcome_node = XY[["outcome_node"]]

  ### Step3) EDNN training ###
  ## Process:
  ## 1. Train EDNN
  ## 2. Calculate the embedding space
  ## To have a stable measure, we repeat the process several times,
  ## which results in several distance vectors in the next step

  embeddingSpaceList = list()
  for (rep in 1:train.rep)
    embeddingSpaceList[[rep]] = EDNN(X ,Y, Xtest = X,
                                     embedding.size, epochs, batch.size, l2reg,
                                     demo = demo, verbose = verbose)

  embeddingSpaceList[["outcome_node"]] = outcome_node
  embeddingSpaceList[["label_node"]] = colnames(X)
  return(embeddingSpaceList)
}


#' Detecting the nodes whose local neighbors change bweteen the two conditions.
#'
#' @param embeddingSpaceList a list obtained by the \code{mnda_embedding_2layer()} function.
#' @param p.adjust.method method for adjusting p-value (including methods on \code{p.adjust.methods}).
#' If set to "none" (default), no adjustment will be performed.
#' @param alpha numeric value of significance level (default: 0.05)
#' @param rank.prc numeric value of the rank percentage threshold (default: 0.1)
#' @param volcano.plot boolean value for generating the Volcano plot (default: TRUE)
#' @param ranksum.sort.plot boolean value for generating the sorted rank sum plot (default: FALSE)
#'
#' @return the highly variable nodes
#' @export
#'
#' @details
#' Calculating the distance of node pairs in the embedding space and check their significance.
#' To find the significantly varying nodes in the 2-layer-network, the distance between
#' the corresponding nodes are calculated along with the null distribution.
#' The null distribution is obtained based on the pairwise distances on null graphs.
#' if in \code{mnda_embedding_2layer} function \code{null.perm=FALSE}, the multiplex network
#' does not have the two randomly permuted graphs, thus the distances between all the nodes will
#' be used for the null distribution.
#'
#' @examples
#' myNet = network_gen(N_nodes = 50, N_var_nodes = 5, N_var_nei = 40, noise_sd = .01)
#' graph_data = myNet[["data_graph"]]
#' embeddingSpaceList = mnda_embedding_2layer(graph.data=graph_data, train.rep=5, walk.rep=5)
#' Nodes = mnda_node_detection_2layer(embeddingSpaceList)
#'
mnda_node_detection_2layer = function(embeddingSpaceList, p.adjust.method = "none",
                                      alpha = 0.05, rank.prc = .1,
                                      volcano.plot = TRUE, ranksum.sort.plot = FALSE){

  node_labels = embeddingSpaceList[["label_node"]]

  outcome_node = embeddingSpaceList[["outcome_node"]]
  outcome = unique(outcome_node)

  assertthat::assert_that(length(outcome) == 4 | length(outcome) == 2)
  if (length(outcome) == 4)
    mode = "null_perm"
  else
    mode = "null_simple"

  N_nodes = nrow(embeddingSpaceList[[1]]) / length(outcome)


  Rep = length(embeddingSpaceList)-2     #the "embeddingSpaceList" consists of two extra elements
  Dist = matrix(0, Rep, N_nodes)
  Dist_null = matrix(0, Rep, N_nodes)
  colnames(Dist) = node_labels

  ### P-value analysis
  ### Calculate significancy p-values of distances in comparison with the null modele
  ### and aggregate p-values across different repeats
  P_value = matrix(0, Rep, N_nodes)
  for (rep in 1:Rep){
    embeddingSpace = embeddingSpaceList[[rep]]
    embeddingSpace_1 = embeddingSpace[outcome_node==outcome[1], ]
    embeddingSpace_2 = embeddingSpace[outcome_node==outcome[2], ]
    for (i in 1:N_nodes)
      Dist[rep,i] = Distance(embeddingSpace_1[i,], embeddingSpace_2[i,], method = "cosine")

    ## Null distribution:
    if (mode == "null_perm"){
      embeddingSpace_1 = embeddingSpace[outcome_node==outcome[3], ]
      embeddingSpace_2 = embeddingSpace[outcome_node==outcome[4], ]
      for (i in 1:N_nodes)
        Dist_null[rep,i] = Distance(embeddingSpace_1[i,], embeddingSpace_2[i,], method = "cosine")
    }else if (mode == "null_simple"){
      embeddingSpace_1 = embeddingSpace[outcome_node==outcome[1], ]
      embeddingSpace_2 = embeddingSpace[outcome_node==outcome[2], ]
      for (i in 1:N_nodes){
        ii = sample(1:N_nodes, 2, replace = FALSE)
        Dist_null[rep,i] = Distance(embeddingSpace_1[ii[1],], embeddingSpace_2[ii[2],], method = "cosine")
      }
    }
    for (i in 1:N_nodes)
      P_value[rep,i] = p_val_rank(x = Dist[rep,i], null.pdf = Dist_null[rep,],
                                  alternative = "greater")
  }

  P_value_aggr = apply(P_value, 2, aggregation::fisher)
  Q_value_aggr = stats::p.adjust(P_value_aggr, method = p.adjust.method)

  significant_nodes_index = which(Q_value_aggr < alpha)
  significant_nodes = node_labels[significant_nodes_index]

  Q_value_aggr_log = -log(Q_value_aggr)
  alpha_log = -log(alpha)

  ### Distance rank sum analysis
  ### Rank node distances and calculate the rank sum across different repeats
  ### Here the "decreasing = FALSE", which means that higher distances correspond to
  ### higher ranks.
  Rank_dist = matrix(0, Rep, N_nodes)
  colnames(Rank_dist) = node_labels
  for (rep in 1:Rep)
    Rank_dist[rep,] = Rank(Dist[rep,], decreasing = FALSE)

  Rank_sum_dist = apply(Rank_dist, 2, sum)

  max_rank_sum_dist = Rep*N_nodes
  rank_threshold = max_rank_sum_dist - rank.prc*max_rank_sum_dist

  high_ranked_nodes_index = which(Rank_sum_dist > rank_threshold)
  high_ranked_nodes = node_labels[high_ranked_nodes_index]

  ### Volcano Plot
  if (volcano.plot){
    plot(Rank_sum_dist, Q_value_aggr_log, pch =  20,
         xlab = "Distance rank sum", ylab = "-log(p.values)")
    graphics::abline(h = alpha_log, col = "red")
    graphics::abline(v = rank_threshold, col = "red")
  }

  if (ranksum.sort.plot){
    colours = rep("grey", N_nodes)
    colours[N_nodes:(N_nodes-length(high_ranked_nodes_index))] = "red"
    plot(sort(Rank_sum_dist), pch = 20, col = colours)
  }

  high_var_nodes_index = intersect(significant_nodes_index, high_ranked_nodes_index)
  high_var_nodes = node_labels[high_var_nodes_index]

  Results = list()
  Results[["node_labels"]] = node_labels
  Results[["p_values"]] = Q_value_aggr
  Results[["rank_sum_dist"]] = Rank_sum_dist

  Results[["significant_nodes_index"]] = significant_nodes_index
  Results[["high_ranked_nodes_index"]] = high_ranked_nodes_index
  Results[["high_var_nodes_index"]] = high_var_nodes_index

  Results[["significant_nodes"]] = significant_nodes
  Results[["high_ranked_nodes"]] = high_ranked_nodes
  Results[["high_var_nodes"]] = high_var_nodes

  return(Results)
}
