
#' Multiplex Network Generation
#'
#' @param N_nodes number of nodes in the graph
#' @param N_var_nodes number of nodes whose neighborhood should change from layer 1 to 2
#' @param N_var_nei number of neighbors that should be changing from layer 1 to 2
#' @param noise_sd the standard deviation of the noise added to the edge weights
#'
#' @return No return value, called to plot subgraphs
#' @export
#'
#' @details In this script we generate random pairs of gene co-expression networks,
#' which are different only in a few (pre-set) number of nodes.
#'
#' @examples
#' myNet = network_gen(N_nodes = 100)
#' graphData = myNet[["data_graph"]]
#' varNodes = myNet[["var_node_set"]]
#'
network_gen = function(N_nodes=100, N_var_nodes=5, N_var_nei=90, noise_sd=.1){

  assertthat::assert_that(N_var_nodes <= N_nodes)
  if (N_var_nei >= N_nodes)
    N_var_nei = N_nodes - 1

  Adj1 = matrix(stats::runif(N_nodes^2), N_nodes)
  Adj2 = Adj1

  var_node_set = sample(N_nodes, N_var_nodes, replace = FALSE)
  for (node in var_node_set){
    var_nei_set = sample((1:N_nodes)[-node], N_var_nei, replace = FALSE)
    new_edge_weights = stats::runif(N_var_nei)
    Adj2[node, var_nei_set] = new_edge_weights
    Adj2[var_nei_set, node] = Adj2[node, var_nei_set]
  }

  noise = abs(matrix(stats::rnorm(N_nodes^2, 0, noise_sd), N_nodes))
  Adj2 = Adj2 + noise

  diag(Adj1) = 0
  diag(Adj2) = 0
  Adj1[lower.tri(Adj1)] = t(Adj1)[lower.tri(Adj1)]
  Adj2[lower.tri(Adj2)] = t(Adj2)[lower.tri(Adj2)]


  EdgeList = data.frame(V1 = as.character(rep(1:N_nodes, times = N_nodes)),
                        V2 = as.character(rep(1:N_nodes, each = N_nodes)))
  EdgeWeights = data.frame(W1 = as.numeric(Adj1),
                           W2 = as.numeric(Adj2))
  data_graph = cbind(EdgeList,EdgeWeights)

  Result = list()
  Result[["data_graph"]] = data_graph
  Result[["var_nodes"]] = var_node_set
  return(Result)
}
