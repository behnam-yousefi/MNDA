
#' Visualization of a subgroup using a circular graph
#'
#' @param graph an igraph object
#' @param node_set the names or indices of the nodes around which the subgroup is plotted.
#' @param labels the labels of the nodes to be indicated. Labels should be a named vector if the \code{node_set} consists of the node names.
#' @param node.importance named numeric vector of the node importance to sort the nodes clockwise.
#' @param n.nodes number of nodes to be displayed. If NULL, all the \code{node_set} and their neighbors are considered.
#' @param node_size size of the nodes in plot (default: 5)
#' @param font_size font size of labels if available (default: 4)
#' @param edge_width numeric value to adjust the thickness of the edges in plot.
#' Two modes are defined: [i] two numbers indicating the min and max (default: c(0.5,4));
#' or [ii] a single number that weights the min/max of original edge weights.
#' @param margin the figure margin (default: 2.5)
#'
#' @return nothing to return
#' @export
#'
#' @details This function plots a sub-graph given by a set of nodes as circular plot.
#' the main inputs to the function are: a graph (as an igraph object) and a set of nodes
#' (e.g. highly variable nodes) around which the subgroup is calculated.
#'
#' @examples
#' data = example_data()
#' subgraph_plot(graph = data[["igraph_example"]], node_set = "a")
#'
subgraph_plot = function(graph, node_set, labels=NULL, node.importance = NULL, n.nodes = NULL,
                         node_size=5, font_size=4, edge_width=c(.5,4), margin=2.5){

  if (is.null(labels)){
    labels = names(igraph::V(graph))
    names(labels) = labels
  }else if (is.null(names(labels))){
    names(labels) = names(igraph::V(graph))
  }

  if (is.null(n.nodes))
    n.nodes = Inf

  # obtain the neighbors of the input node_set
  node_set_neigh = c()
  for (n in node_set)
    node_set_neigh = c(node_set_neigh, names(igraph::neighbors(graph, n)))
  Nodes = unique(c(node_set, node_set_neigh))

  if (!is.null(node.importance))
    Nodes = names(sort(node.importance[Nodes], decreasing = TRUE))[1:min(n.nodes, length(Nodes))]

  N_node_set = length(node_set)
  # N_node_neigh = length(node_set_neigh)
  N_nodes = length(Nodes)

  # to obtain the subgraph:
  # [i] calculate the the adjacency matrix of the subgraph
  # [ii] then build the subgraph object
  Adj_mat = as.matrix(igraph::as_adj(graph,  attr = "weight"))
  Adj_mat_sub = Adj_mat[Nodes, Nodes]
  graph_sub = igraph::graph_from_adjacency_matrix(Adj_mat_sub, weighted = TRUE, mode = "undirected")
  W = igraph::E(graph_sub)$weight
  labels = labels[Nodes]

  # set colors and fonts
  EdgeColor = ifelse(W>0, "hotpink", "seagreen2")
  NodeColor = rep("slategray2", length(Nodes))
  NodeColor[1:N_node_set] = "slategrey"
  TxtColor = rep("slategrey", length(Nodes))
  TxtColor[1:N_node_set] = "black"
  FontFace = rep("plain", length(Nodes))
  FontFace[1:N_node_set] = "plain"
  FontSize = rep(font_size, length(Nodes))
  FontSize[1:N_node_set] = font_size
  if (length(edge_width)>1){
    edgeWidth = edge_width
  }else{
    edgeWidth = c(min(abs(W)), max(abs(W))) * edge_width
  }
  # parameters used to account for the circulatory text problem
  a = floor(180 / (360 / N_nodes))
  b = N_nodes - a

  x = 1
  y = 1
  # final plot
  # requireNamespace("ggplot2")
  # requireNamespace("ggraph")
  ggraph::ggraph(graph_sub, layout = 'linear', circular = TRUE) +
    ggraph::geom_edge_arc(ggplot2::aes(color = EdgeColor, edge_width = abs(W)), edge_alpha = .7) +
    ggraph::geom_node_point(ggplot2::aes(x = x, y = y), size = node_size, fill = NodeColor,
                    shape = 21, colour = "black", stroke = .4) +
    ggraph::geom_node_text(ggplot2::aes(x = x*1.12, y = y*1.12,
                       angle = -((-ggraph::node_angle(x, y)+90)%%180)+90, hjust = c(rep(0,a), rep(1,b))),
                   label = labels,
                   fontface = FontFace,
                   size = FontSize,
                   color = TxtColor) +
    ggraph::scale_edge_width_continuous(range = edgeWidth) +
    ggplot2::scale_size_identity() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position="none") +
    ggplot2::expand_limits(x = c(-margin, margin), y = c(-margin, margin))
}



#' Visualization of a difference subgroup using a circular graph
#'
#' @param mnda.graph mnda.graph data
#' @param node.importance named numeric vector of the node importance to sort the nodes clockwise.
#' @param n.var.nodes number of variable nodes to show
#' @param n.neigh number of neighboring  nodes to show
#' @param diff.threshold edge threshold
#' @param edge.width numeric value to adjust the thickness of the edges in plot.
#' Two modes are defined: [i] two numbers indicating the min and max (default: c(0.5,4));
#' or [ii] a single number that weights the min/max of original edge weights.
#'
#' @return nothing to return
#' @export
#'
#'#' @details This function plots a difference sub-graph as circular plot.
#' the main inputs to the function are: a graph multiplex (as an mnda.graph) and a vector of node importances.
#'
#' @examples
#' myNet = network_gen(N_nodes = 100, N_var_nodes = 5, N_var_nei = 90, noise_sd = .01)
#' graph_data = myNet[["data_graph"]]
#' node_importance_dummy = 1:100
#' names(node_importance_dummy) = 1:100
#' subgraph_difference_plot(graph_data, node.importance = node_importance_dummy)
#'
subgraph_difference_plot = function(mnda.graph, node.importance,
                                    n.var.nodes = 5, n.neigh = 10,
                                    diff.threshold=0, edge.width=c(.5,4)){

  var_nodes = sort(node.importance, decreasing = TRUE)[1:n.var.nodes]
  var_nodes = names(var_nodes)

  graph_to_plot = cbind(mnda.graph[,1:2],
                        W = as.numeric(mnda.graph[,3]) - as.numeric(mnda.graph[,4]))
  G = mnda::as.igraph(graph_to_plot, diff.threshold)

  # hist(graph_to_plot$W)
  # hist(E(G)$weight)

  subgraph_plot(G, var_nodes, node.importance = node.importance,
                n.nodes = (n.var.nodes + n.neigh), edge_width = edge.width)

}
