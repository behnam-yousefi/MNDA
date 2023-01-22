## Pipeline for multiplex network differential analysis
## The aim of the current functions:
## Perform MNDA on multiple networks (e.g. individual specific networks - ISNs)
## with two outcomes and test whether the distances are significant.
## Networks are considered to be paired with respect to the outcome.
## Example of the usage is provided in ???

# Null hypothesis (H0): The distributions of pairwise distances across the two groups of ISNs
# are not significantly different.
## The test is performed by either t.test or wilcox.test.

#' Calculate the embedding space for a multiplex network
#'
#' @param graph.data dataframe of the graph data containing edge list and edge weights.
#' column 1 and 2 consisting of the edge list (undirected).
#' column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
#' @param outcome a vector of outcomes for each network.
#' @param indv.index the index of individual networks.
#' @param edge.threshold numeric value to set edge weights below the threshold to zero (default: 0). the greater edge weights do not change.
#' @param train.rep numeric value to set the number of EDNN training repeats (default: 50).
#' @param embedding.size the dimension of embedding space, equal to the number of the bottleneck hidden nodes (default: 5).
#' @param epochs maximum number of pocks. An early stopping callback with a patience of 5 has been set inside the function (default = 10).
#' @param batch.size batch size for learning (default = 5).
#' @param l2reg the coefficient of L2 regularization for the input layer (default = 0).
#' @param walk.rep number of repeats for the random walk (default: 100).
#' @param n.steps number of the random walk steps (default: 5).
#' @param random.walk boolean value to enable the random walk algorithm (default: TRUE).
#' @param demo a boolean vector to indicate this is a demo example or not
#' @param verbose if \emph{TRUE} a progress bar is shown.
#'
#' @return a list of embedding spaces for each node.
#' @export
#'
#' @examples
#' myNet = network_gen(N_nodes = 50, N_var_nodes = 5, N_var_nei = 40, noise_sd = .01)
#' graph_data = myNet[["data_graph"]]
#' embeddingSpaceList = mnda_embedding(graph.data=graph_data, outcome=c(1,2),
#' train.rep=2, random.walk=FALSE)
#'
mnda_embedding = function(graph.data, outcome, indv.index = NULL,
                          edge.threshold=0, train.rep=50,
                          embedding.size=2, epochs=10, batch.size=5, l2reg=0,
                          walk.rep=100, n.steps=5, random.walk=TRUE, demo = TRUE, verbose=TRUE){

  assertthat::assert_that(train.rep >= 2)

  ### Step1) Read Graph Data ###
  ## The graph.data format is a data.frame:
  ## column 1 and 2 consisting of the edge list (undirected)
  ## column 3 and 4 consisting the edge weights corresponding to each graph, respectively.
  EdgeList = graph.data[,1:2]
  EdgeWeights = graph.data[,-c(1,2)]

  N_network = ncol(EdgeWeights)

  ### Step2) Preparing the input and output of the EDNN for all the graphs ###
  XY = ednn_IOprepare(edge.list = EdgeList, edge.weight = EdgeWeights, indv.index = indv.index,
                      walk.rep = walk.rep, n.steps = n.steps, random.walk = random.walk,
                      outcome = outcome, edge.threshold = edge.threshold, verbose = verbose)
  X = XY[["X"]]
  Y = XY[["Y"]]
  outcome_node = XY[["outcome_node"]]
  individual_node = XY[["individual_node"]]

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
  embeddingSpaceList[["individual_node"]] = individual_node
  embeddingSpaceList[["label_node"]] = colnames(X)
  return(embeddingSpaceList)
}


#' Detecting the nodes whose local neighbors change between the two conditions for ISNs.
#'
#' @param embeddingSpaceList a list obtained by the \code{mnda_embedding_2layer()} function.

#' @return the distances for each repeat
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
#' embeddingSpaceList = mnda_embedding(graph.data=graph_data, outcome=c(1,2),
#' indv.index=c(1,1), train.rep=2, random.walk=FALSE)
#' Dist = mnda_node_distance(embeddingSpaceList)
#'
mnda_node_distance = function(embeddingSpaceList){

  node_labels = embeddingSpaceList[["label_node"]]
  N_nodes = length(node_labels)

  outcome_node = embeddingSpaceList[["outcome_node"]]
  outcome = unique(outcome_node)
  assertthat::assert_that(length(outcome) == 2)

  individual_node = embeddingSpaceList[["individual_node"]]
  individual = unique(individual_node)
  N_indv = length(individual)

  Rep = length(embeddingSpaceList)-3     #the "embeddingSpaceList" consists of three extra elements
  Dist_list = list()

  ### P-value analysis
  ### Calculate significancy p-values of distances in comparison with the null modele
  ### and aggregate p-values across different repeats
  for (rep in 1:Rep){
    Dist = matrix(0, N_indv, N_nodes)
    rownames(Dist) = individual
    colnames(Dist) = node_labels

    embeddingSpace = embeddingSpaceList[[rep]]

    for (indv in 1:N_indv){
      embeddingSpace_1 = embeddingSpace[individual_node==individual[indv] & outcome_node==outcome[1], ]
      embeddingSpace_2 = embeddingSpace[individual_node==individual[indv] & outcome_node==outcome[2], ]
      for (i in 1:N_nodes)
        Dist[indv,i] = Distance(embeddingSpace_1[i,], embeddingSpace_2[i,], method = "cosine")
    }
    Dist_list[[rep]] = Dist
  }

  return(Dist_list)
}

#' Test the embedding distances of local neighbors change between the two conditions for ISNs.
#'
#' @param Distance a distance list obtained by the \code{mnda_node_distance()} function.
#' @param y vector with the length equal to the number of individuals.
#' @param stat.test statistical test used to detect the nodes \code{c("t.test","wilcox.test")} (default: wilcox.test)
#' @param p.adjust.method method for adjusting p-value (including methods on \code{p.adjust.methods}).
#' If set to "none" (default), no adjustment will be performed.
#'
#' @return The adjusted pvalues for each node.
#' @export
#'
#' @details
#' The adjusted p-values for each node is calculated based on their distance variation between
#' the two conditions.
#'
#' @examples
#' ISN1 = network_gen(N_nodes = 50, N_var_nodes = 5, N_var_nei = 40, noise_sd = .01)
#' ISN2 = network_gen(N_nodes = 50, N_var_nodes = 5, N_var_nei = 40, noise_sd = .01)
#' graph_data = cbind(ISN1[["data_graph"]], ISN1[["data_graph"]][,3:4])
#' embeddingSpaceList = mnda_embedding(graph.data=graph_data, outcome=c(1,2,1,2),
#' indv.index=c(1,1,2,2), train.rep=2, random.walk=FALSE)
#' Dist = mnda_node_distance(embeddingSpaceList)
#' Result = mnda_distance_test_isn(Dist, y = c(1,2))
#'
mnda_distance_test_isn  = function(Distance, y,
                               stat.test = "wilcox.test", p.adjust.method = "none"){

  Rep = length(Distance)
  N_nodes = ncol(Distance[[1]])

  y_unique = unique(y)
  assertthat::assert_that(length(y_unique)==2)

  P_value = matrix(0, Rep, N_nodes)

  for (rep in 1:Rep){
    Dist = Distance[[rep]]

    for (i in 1:N_nodes){
      x1 = Dist[y == y_unique[1],i]
      x2 = Dist[y == y_unique[2],i]

      if (stat.test == "wilcox.test")
        p.val = stats::wilcox.test(x1, x2)$p.val
      else if (stat.test == "t.test")
        p.val = stats::t.test(x1, x2)$p.val
      else
        message("error: Not a valid stat.test value. Possible values: wilcox.test, t.test")

      P_value[rep, i] = p.val
    }
  }

  P_value_aggr = apply(P_value, 2, aggregation::fisher)
  Q_value_aggr = stats::p.adjust(P_value_aggr, method = p.adjust.method)

  names(Q_value_aggr) = colnames(Distance[[1]])
  return(Q_value_aggr)
}
