
#' Visualization of a subgroup using a circular graph
#'
#' @param graph an igraph object
#' @param node_set the names or indices of the nodes around which the subgroup is plotted.
#' @param labels the labels of the nodes to be indicated. Labels should be a named vector if the \code{node_set} consists of the node names.
#' @param node_size size of the nodes in plot (default: 5)
#' @param font_size font size of labels if available (default: 4)
#' @param max_edge_width maximum thickness of the edges in plot (default: 4)
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
#'subgraph_plot(graph = igraph_example, node_set = "a")
#'
subgraph_plot = function(graph, node_set, labels=NULL,
                         node_size=5, font_size=4, max_edge_width=4, margin=2.5){

  if (is.null(labels)){
    labels = names(igraph::V(graph))
    names(labels) = labels
  }else if (is.null(names(labels))){
    names(labels) = names(igraph::V(graph))
  }

  # obtain the neighbors of the input node_set
  node_set_neigh = c()
  for (n in node_set)
    node_set_neigh = c(node_set_neigh, names(igraph::neighbors(graph, n)))
  Nodes = unique(c(node_set, node_set_neigh))

  N_node_set = length(node_set)
  N_node_neigh = length(node_set_neigh)
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

  # parameters used to account for the circulatory text problem
  a = floor(180 / (360 / N_nodes))
  b = N_nodes - a

  # final plot
  requireNamespace(ggplot2)
  requireNamespace(ggraph)
  ggraph(graph_sub, layout = 'linear', circular = TRUE) +
    geom_edge_arc(aes(color = EdgeColor, edge_width = abs(W)^1), edge_alpha = .7) +
    geom_node_point(aes(x = x, y = y), size = node_size, fill = NodeColor,
                    shape = 21, colour = "black", stroke = .4) +
    geom_node_text(aes(x = x*1.12, y = y*1.12,
                       angle = -((-node_angle(x, y)+90)%%180)+90, hjust = c(rep(0,a), rep(1,b))),
                   label = labels,
                   fontface = FontFace,
                   size = FontSize,
                   color = TxtColor) +
    scale_size_identity() +
    scale_edge_width_continuous(range = c(.0001, max_edge_width)) +
    theme_void() +
    theme(legend.position="none") +
    expand_limits(x = c(-margin, margin), y = c(-margin, margin))


}

