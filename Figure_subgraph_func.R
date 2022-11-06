## Visualization of a subgroup using a circular graph
## The input is a set of nodes (e.g. highly variable nodes),
## around which the subgroup is calculated.

library(igraph)
library(ggraph)

subgraph_plot = function(graph, node_set, labels=NULL, margin=2.5){
  # graph: An igraph object
  # node_set: The names or indices of the nodes around which the subgroup is plotted.
  # labels: The labels of the nodes to be indicated. labels should be named if the node_set 
  # consists of the node names.
  # margin: the figure margin.
  
  if (is.null(labels)){
    labels = names(V(graph))
    names(labels) = labels
  }else if (is.null(names(labels))){
    names(labels) = names(V(graph))
  }
  
  # obtain the neighbors of the input node_set
  node_set_neigh = c()
  for (n in node_set)
    node_set_neigh = c(node_set_neigh, names(neighbors(graph, n)))
  Nodes = unique(c(node_set, node_set_neigh))
  
  N_node_set = length(node_set)
  N_node_neigh = length(node_set_neigh)
  N_nodes = length(Nodes)
  
  # to obtain the subgraph: 
  # [i] calculate the the adjacency matrix of the subgraph 
  # [ii] then build the subgraph object
  Adj_mat = as.matrix(as_adj(graph,  attr = "weight"))
  Adj_mat_sub = Adj_mat[Nodes, Nodes]
  graph_sub = graph_from_adjacency_matrix(Adj_mat_sub, weighted = TRUE, mode = "undirected")
  W = E(graph_sub)$weight
  labels = labels[Nodes]
  
  # set colors and fonts
  EdgeColor = ifelse(W>0, "hotpink", "seagreen2")
  NodeColor = rep("slategray2", length(Nodes))
  NodeColor[1:N_node_set] = "slategrey"
  TxtColor = rep("slategrey", length(Nodes))
  TxtColor[1:N_node_set] = "black"
  FontFace = rep("plain", length(Nodes))
  FontFace[1:N_node_set] = "plain"
  FontSize = rep(4, length(Nodes))
  FontSize[1:N_node_set] = 4
  
  # parameters used to account for the circulatory text problem
  a = floor(180 / (360 / N_nodes))
  b = N_nodes - a
  
  # final plot
  ggraph(graph_sub, layout = 'linear', circular = TRUE) +
    geom_edge_arc(aes(color = EdgeColor, edge_width = abs(W)^1), edge_alpha = .7) +
    geom_node_point(aes(x = x, y = y), size = 5, fill = NodeColor, 
                    shape = 21, colour = "black", stroke = .4) +
    geom_node_text(aes(x = x*1.12, y = y*1.12, 
                  angle = -((-node_angle(x, y)+90)%%180)+90, hjust = c(rep(0,a), rep(1,b))), 
                  label = labels,
                  fontface = FontFace,
                  size = FontSize, 
                  color = TxtColor) +
    scale_size_identity() +
    scale_edge_width_continuous(range = c(.0001, 3)) + 
    theme_void() +
    theme(legend.position="none") + 
    expand_limits(x = c(-margin, margin), y = c(-margin, margin))
  
  
}
  
