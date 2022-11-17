## Multiplex Graph Representation Learning and
## multiplex network differential analysis
## By: Behnam Yousefi

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


#' Encoder decoder neural network (EDNN) function
#'
#' @param X concatenated adjacency matrices for different layers containing the nodes in training phase
#' @param Y concatenated random walk probability matrices for different layers containing the nodes in training phase
#' @param Xtest concatenated adjacency matrices for different layers containing the nodes in test phase. Can be = \emph{X} for transductive inference.
#' @param embedding_size the dimension of embedding space, equal to the number of the bottleneck hidden nodes.
#' @param epochs maximum number of pocks. An early stopping callback with a patience of 5 has been set inside the function (default = 10).
#' @param batch_size batch size for learning (default = 5).
#' @param l2reg the coefficient of L2 regularization for the input layer (default = 0).
#'
#' @return The embedding space for Xtest.
#' @export
#'
#' @examples
#' embeddingSpace = EDNN(X = Adj, Y = Prob, Xtest = Adj)
#'
EDNN = function(X, Y, Xtest, embedding_size = 2, epochs = 10, batch_size = 5, l2reg = 0){

  Nnode = ncol(X)
  inputSize = ncol(X)
  outputSize = ncol(Y)

  # Define Encoder
  enc_input = layer_input(shape = inputSize)
  enc_output = enc_input %>%
    layer_dense(units=embedding_size, activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg))

  encoder = keras_model(enc_input, enc_output)

  # Define decoder
  dec_input = layer_input(shape = embedding_size)
  dec_output = dec_input %>%
    layer_dense(units = outputSize, activation = "sigmoid")

  decoder = keras_model(dec_input, dec_output)

  # Define Auto-Encoder
  aen_input = layer_input(shape = inputSize)
  aen_output = aen_input %>%
    encoder() %>%
    decoder()

  autoencoder = keras_model(aen_input, aen_output)

  # Training configuration
  autoencoder %>% compile(loss = "mse", optimizer = 'adam')
  checkpoint <- callback_model_checkpoint(filepath = "My_model_temp", save_best_only = TRUE, verbose = 0)
  early_stopping <- callback_early_stopping(patience = 5)

  # Fit the model and save in history
  history <- autoencoder %>% fit(X, Y, validation_data = list(X, Y), loss = "mse",
                                 epochs = epochs, batch_size = batch_size, callbacks = list(checkpoint, early_stopping))


  # Final embeding
  embeddingSpace = encoder %>% predict(Xtest)
  return(embeddingSpace)
}


### function to calculate distance between two vectors
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

#' Preparing the input and output of the EDNN for a multiplex graph
#'
#' @param edge.list edge list as a dataframe with two columns.
#' @param edge.weight edge weights as a dataframe. Each column corresponds to a graph. By default, the \code{colnames} are considered as outcomes unless indicated in \code{outcome} argument.
#' @param outcome clinical outcomes for each graph. If not mentioned, the \code{colnames(edge.weight)} are considered by default.
#' @param edge.threshold numeric value to set edge weights below the threshold to zero (default: 0). the greater edge weights do not change.
#' @param walk.rep number of repeats for the random walk (default: 100).
#' @param n.steps number of the random walk steps (default: 5).
#'
#' @return the input and output required to train the EDNN
#' @export
#'
#' @examples
#' XY = ednn_IOprepare(edge.list, edge.weight)
#' embedding = EDNN(X = XY["X"] ,Y = XY["Y"], Xtest = XY["X"])
#'
ednn_IOprepare = function(edge.list, edge.weight, outcome=NULL, edge.threshold=0,
                          walk.rep=100, n.steps=5){

  if (is.null(outcome))
    outcome = colnames(edge.weight)

  graph = graph(t(edge.list), directed = FALSE)

  N_nodes = length(V(graph))
  N_graph = ncol(edge.weight)

  X = c()
  Y = c()
  outcome_node = c()

  for (i in 1:N_graph){
    ## Set layer-specific weights for each graph and
    ## [optional] perform a thresholding if needed
    W_i = as.numeric(edge.weight[,i])
    graph_i = simplify(set.edge.attribute(graph, "weight", index=E(graph), W_i))
    graph_i = delete_edges(graph_i, E(graph_i)[E(graph_i)$weight < edge.threshold])

    ## Step 1) Input: Adjacency matrix calculation
    Adj_i = as.matrix(as_adj(graph_i,  attr = "weight"))

    ## Step 2) Output: Perform the fixed-length random walk
    ## and calculating the node visit probabilities.
    # Two options exist:
    # 1.repetitive simple random walks
    # 2.repetitive weighted random walks (specific to weighted graphs)
    RW = RepRandomWalk (graph_i, Nrep = walk.rep, Nstep = n.steps, weighted_walk = TRUE)

    ## Step 3) Make it multilayer for EDNN
    X = rbind(X, Adj_i)
    Y = rbind(Y, RW$Probabilities)

    outcome_node = c(outcome_node, rep(outcome[i], N_nodes))
  }
  X = X / (apply(X, 1, sum) + .000000001)

  return(list(X,Y))
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


