## Multiplex Graph Representation Learning and
## multiplex network differential analysis

#' Function to calculate distance between two vectors
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param method distance calculation method: cosine (default), dot.prod, euclidian, manhattan, chebyshev, coassociation
#'
#' @return the distance value
#' @export
#'
#' @examples
#' x = c(1,2,3)
#' y = c(6,4,6)
#' Distance(x,y)
#'
Distance = function(x, y, method = "cosine"){
  dist = switch(method,
                "cosine" = 1- ((x %*% y) / sqrt((sum(x^2))*(sum(y^2)) + .000000001)),
                "dot.prod" = 1/(x %*% y),
                "euclidian" = sum((x - y)^2),
                "manhattan" = sum(abs(x - y)),
                "chebyshev" = max(abs(x-y)),
                # "coassociation" = 1 - coassociation_sim(x,y)
  )
  return(dist)
}


#' Ranking a vector
#'
#' @param x a numeric vector
#' @param decreasing logical. Should the sort order be increasing or decreasing? (defualt: FALSE)
#'
#' @return the rank of the vector elements
#' @export
#'
#' @details
#' hint: What is the difference between Order and Rank\cr
#' Order: [the index of the greatest number, ..., the index of the smallest number]\cr
#' Rank: [the rank of the 1st number, ..., the rank of the last number]\cr
#' In Rank, the order of the numbers remains constant so can be used for ranksum.\cr
#' ex) \cr
#'     > a = c(10, 20, 50, 30, 40)\cr
#'     > order(a)\cr
#'     [1] 1 2 4 5 3]]\cr
#'     > Rank(a)\cr
#'     [1] 1 2 5 3 4
#'
#' @examples
#' a = c(10, 20, 50, 30, 40)
#' Rank(a)
#'
Rank = function(x, decreasing = FALSE) {

  Order = order(x, decreasing = decreasing)

  Rank = rep(0,length(Order))
  for(i in 1:length(Order))
    Rank[Order[i]] = i

  return(Rank)
}

#' Convert adjacency matrix to mnda graph data
#'
#' @param adj.list list of adjacency matrices with matching nodes
#' @param outcome graph outcomes or graph labels. If NULL, \code{outcome = 1:N_graphs}.
#'
#' @return mnda.graph data
#' @export
#'
#' @examples
#' data = example_data()
#' adj.list = list(data[["adj_mat_example"]], data[["adj_mat_example"]])
#' graph.data = as.mnda.graph(adj.list)
#'
as.mnda.graph = function(adj.list, outcome = NULL){

  N_graphs = length(adj.list)
  N_nodes = ncol(adj.list[[1]])

  if (is.null(outcome))
    outcome = 1:N_graphs

  node_labels = rownames(adj.list[[1]])
  if (is.null(node_labels))
    node_labels = 1:N_nodes

  ## Set edge list
  EdgeList = data.frame(V1 = rep(node_labels, times = N_nodes),
                        V2 = rep(node_labels, each = N_nodes))

  ## Set edge weight
  EdgeWeights = c()
  for (i in 1:N_graphs)
    EdgeWeights = cbind(EdgeWeights, as.numeric(adj.list[[i]]))
  colnames(EdgeWeights) = outcome

  ## set mnda graph data
  data_graph = cbind(EdgeList,EdgeWeights)

  return(data_graph)
}

#' Convert mnda graph data to igraph
#'
#' @param mnda.graph mnda graph data
#' @param edge.threshold numeric
#'
#' @return igraph object
#' @export
#'
#' @examples
#' data = example_data()
#' graph = as.igraph(mnda.graph = data[["mnda_graph_example"]])
#'
as.igraph = function(mnda.graph, edge.threshold=0){

  EdgeList = mnda.graph[,1:2]
  EdgeWeights = as.numeric(mnda.graph[,3])

  EdgeList = t(EdgeList)
  graph = igraph::graph(EdgeList, directed = FALSE)

  graph = igraph::simplify(igraph::set.edge.attribute(graph, "weight", index=igraph::E(graph), EdgeWeights/2))   # EdgeWeights/2 is used due to the use of simplify() function
  graph = igraph::delete_edges(graph, igraph::E(graph)[abs(igraph::E(graph)$weight) < edge.threshold])

  return(graph)

}
