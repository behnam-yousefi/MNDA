% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{as_mnda_graph}
\alias{as_mnda_graph}
\title{Convert adjacency matrix to mnda graph data}
\usage{
as_mnda_graph(adj.list, outcome = NULL)
}
\arguments{
\item{adj.list}{list of adjacency matrices with matching nodes}

\item{outcome}{graph outcomes or graph labels. If NULL, \code{outcome = 1:N_graphs}.}
}
\value{
mnda.graph data
}
\description{
Convert adjacency matrix to mnda graph data
}
\examples{
data = example_data()
adj.list = list(data[["adj_mat_example"]], data[["adj_mat_example"]])
graph.data = as_mnda_graph(adj.list)

}
