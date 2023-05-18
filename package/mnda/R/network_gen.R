
#' Multiplex Network Generation
#'
#' @param n.nodes number of nodes in the graph
#' @param n.var.nodes number of nodes whose neighborhood should change from layer 1 to 2
#' @param n.var.nei number of neighbors that should be changing from layer 1 to 2
#' @param noise.sd the standard deviation of the noise added to the edge weights
#'
#' @return No return value, called to plot subgraphs
#' @export
#'
#' @details In this script we generate random pairs of gene co-expression networks,
#' which are different only in a few (pre-set) number of nodes.
#'
#' @examples
#' myNet = network_gen(n.nodes = 100)
#' graphData = myNet[["data_graph"]]
#' varNodes = myNet[["var_node_set"]]
#'
network_gen = function(n.nodes=100, n.var.nodes=5, n.var.nei=90, noise.sd=.1){

  assertthat::assert_that(n.var.nodes <= n.nodes)
  if (n.var.nei >= n.nodes)
    n.var.nei = n.nodes - 1

  Adj1 = matrix(stats::runif(n.nodes^2), n.nodes)
  Adj2 = Adj1

  var_node_set = sample(n.nodes, n.var.nodes, replace = FALSE)
  for (node in var_node_set){
    var_nei_set = sample((1:n.nodes)[-node], n.var.nei, replace = FALSE)
    new_edge_weights = stats::runif(n.var.nei)
    Adj2[node, var_nei_set] = new_edge_weights
    Adj2[var_nei_set, node] = Adj2[node, var_nei_set]
  }

  noise = abs(matrix(stats::rnorm(n.nodes^2, 0, noise.sd), n.nodes))
  Adj2 = Adj2 + noise

  diag(Adj1) = 0
  diag(Adj2) = 0
  Adj1[lower.tri(Adj1)] = t(Adj1)[lower.tri(Adj1)]
  Adj2[lower.tri(Adj2)] = t(Adj2)[lower.tri(Adj2)]


  EdgeList = data.frame(V1 = as.character(rep(1:n.nodes, times = n.nodes)),
                        V2 = as.character(rep(1:n.nodes, each = n.nodes)))
  EdgeWeights = data.frame(W1 = as.numeric(Adj1),
                           W2 = as.numeric(Adj2))
  data_graph = cbind(EdgeList,EdgeWeights)

  Result = list()
  Result[["data_graph"]] = data_graph
  Result[["var_nodes"]] = var_node_set
  return(Result)
}
