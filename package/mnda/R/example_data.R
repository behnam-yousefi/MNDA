#' Example Data
#'
#' @return example data as a list: "adj_mat_example", "igraph_example", "mnda_graph_example"
#' @export
#'
#' @examples
#' data = example_data()

example_data = function(){

  adj_mat_example = rbind(c(0,1,2),c(1,0,3),c(2,3,0))

  igraph_example = igraph::graph_from_adjacency_matrix(adj_mat_example,
                                                       weighted = TRUE,
                                                       mode = "undirected")
  igraph::V(igraph_example)$name = c("a", "b", "c")
  igraph_example = igraph::simplify(igraph::set.edge.attribute(igraph_example, "weight",
                                                               index=igraph::E(igraph_example), c(1,2,3)))


  mnda_graph_example = data.frame(V1 = c(1,1,2),
                                  V2 = c(2,3,3),
                                  W = c(1,2,3))

  data = list()
  data[["adj_mat_example"]] = adj_mat_example
  data[["igraph_example"]] = igraph_example
  data[["mnda_graph_example"]] = mnda_graph_example
  return(data)
}
