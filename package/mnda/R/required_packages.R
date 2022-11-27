# ## Required R packages
#
# library(igraph)
#
# library(MASS)
# library(ggraph)
#
# library(aggregation)
#
# require(Matrix)
# library(keras)
#

# Example Data

Adj_mat_example = rbind(c(0,1,2),c(1,0,3),c(2,3,0))

igraph_example = igraph::graph_from_adjacency_matrix(Adj_mat_example,
                                                     weighted = TRUE,
                                                     mode = "undirected")
igraph::V(igraph_example)$name = c("a", "b", "c")
igraph_example = igraph::simplify(igraph::set.edge.attribute(igraph_example, "weight",
                                                             index=igraph::E(igraph_example), c(1,2,3)))


mnda_graph_example = data.frame(V1 = c(1,1,2),
                                V2 = c(2,3,3),
                                W = c(1,2,3))

