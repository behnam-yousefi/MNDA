#' Weighted Random Walk algorithm
#'
#' @param graph an igraph object
#' @param startNode the starting node (i.e. a node name or a node index)
#' @param maxStep maximum number steps (default:5)
#' @param node_names ???
#'
#' @return The set of nodes passed by the random walker.
#' @export
#'
#' @examples
#' nodePath = WeightdRandomWalk(graph, 1)
#'
WeightdRandomWalk = function(graph, startNode, maxStep = 5, node_names = FALSE){

  A = as.matrix(as_adj(graph, attr = "weight"))
  N_node = ncol(A)
  outputNodes = rep(NA, maxStep)
  nextStep = startNode

  for (i in 1:maxStep) {
    w = A[nextStep,]

    if (sum(w)==0){
      outputNodes = rep(nextStep,maxStep)
      break

    }else{
      w = w/sum(w)
      nextStep = sample(N_node, size = 1, prob = w)
      if (node_names)
        nextStep = rownames(A)[nextStep]
      outputNodes[i] = nextStep
    }
  }
  return(outputNodes)
}


#' Repetitive Fixed-length (weighted) random walk algorithm
#'
#' @param graph an igraph object
#' @param Nrep number of repeats (default:100)
#' @param Nstep maximum number steps (default:5)
#' @param weighted_walk choose the \emph{weighted walk} algorithm if \emph{TRUE} and \emph{simple random} walk if \emph{FALSE}. (default: \emph{TRUE})
#'
#' @return
#' Steps (S): The total number of times a node is visited starting from the corresponding node in the row.
#' Probabilities (P): The node visit probabilities starting from the corresponding node in the row.
#' @export
#'
#' @examples
#' RW = RepRandomWalk(graph)
#' Steps = RW[["Steps"]]
#' Probabilities = RW[["Probabilities"]]
#'
RepRandomWalk = function(graph, Nrep = 100, Nstep = 5, weighted_walk = TRUE){

  N = length(V(graph))
  S = matrix(0, N, N)

  i = 0
  for (node in V(graph)){
    i = i+1

    Neighbors = c(node, neighbors(graph,node))
    NeighborsOneHot = rep(0,N)
    NeighborsOneHot[Neighbors] = 1

    WalkRep = matrix(0, Nrep, N)
    for (rep in 1:Nrep){

      if (weighted_walk)
        # use WeightdRandomWalk() written in the current script.
        RW_steps = WeightdRandomWalk(graph, startNode = node, maxStep = Nstep)
      else
        # use igraph::random_walk
        RW_steps = random_walk(graph, start = node, steps = Nstep)

      WalkRep[rep,RW_steps] = 1
    }

    S[i,] = apply(WalkRep, 2, sum)
  }

  P = S / Nrep
  P[P>1] = 1
  diag(P) = 1

  colnames(S) = names(V(graph))
  colnames(P) = names(V(graph))

  Result = list()
  Result[["Steps"]] = S
  Result[["Probabilities"]] = P

  return(Result)
}