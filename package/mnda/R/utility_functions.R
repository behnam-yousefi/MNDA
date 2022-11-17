## Multiplex Graph Representation Learning and
## multiplex network differential analysis
## By: Behnam Yousefi


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
                "coassociation" = 1 - coassociation_sim(x,y)
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


