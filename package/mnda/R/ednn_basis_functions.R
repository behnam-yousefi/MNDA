#' Encoder decoder neural network (EDNN) function
#' @param X concatenated adjacency matrices for different layers containing the nodes in training phase
#' @param Y concatenated random walk probability matrices for different layers containing the nodes in training phase
#' @param Xtest concatenated adjacency matrices for different layers containing the nodes in test phase. Can be = \emph{X} for transductive inference.
#' @param embedding_size the dimension of embedding space, equal to the number of the bottleneck hidden nodes.
#' @param epochs maximum number of pocks. An early stopping callback with a patience of 5 has been set inside the function (default = 10).
#' @param batch_size batch size for learning (default = 5).
#' @param l2reg the coefficient of L2 regularization for the input layer (default = 0).
#' @param demo a boolean vector to indicate this is a demo example or not
#' @param verbose if \emph{TRUE} a progress bar is shown.
#'
#' @return The embedding space for Xtest.
#' @export
#'
#' @examples
#' myNet = network_gen(N_nodes = 50)
#' graphData = myNet[["data_graph"]]
#' edge.list = graphData[,1:2]
#' edge.weight = graphData[,3:4]
#' XY = ednn_IOprepare(edge.list, edge.weight)
#' X = XY[["X"]]
#' Y = XY[["Y"]]
#' embeddingSpace = EDNN(X = X, Y = Y, Xtest = X)
#'
EDNN = function(X, Y, Xtest, embedding_size = 2, epochs = 10, batch_size = 5, l2reg = 0, demo = TRUE, verbose = TRUE){

  Nnode = ncol(X)
  inputSize = ncol(X)
  outputSize = ncol(Y)

  if (demo){
    return(Y[,1:embedding_size])
  }else{

  # Define Encoder
  enc_input = keras::layer_input(shape = inputSize)
  enc_output = # enc_input %>%
    keras::layer_dense(enc_input, units=embedding_size, activation = "relu",
                kernel_regularizer = keras::regularizer_l2(l2reg))

  encoder = keras::keras_model(enc_input, enc_output)

  # Define decoder
  dec_input = keras::layer_input(shape = embedding_size)
  dec_output = # dec_input %>%
    keras::layer_dense(dec_input, units = outputSize, activation = "sigmoid")

  decoder = keras::keras_model(dec_input, dec_output)

  # Define Auto-Encoder
  aen_input = keras::layer_input(shape = inputSize)
  aen_output = # aen_input %>%
    # encoder() %>%
    decoder(encoder(aen_input))

  autoencoder = keras::keras_model(aen_input, aen_output)

  # Training configuration
  # autoencoder %>% keras::compile(loss = "mse", optimizer = 'adam')
  keras::compile(autoencoder, loss = "mse", optimizer = 'adam')
  checkpoint <- keras::callback_model_checkpoint(filepath = "My_model_temp", save_best_only = TRUE, verbose = 0)
  early_stopping <- keras::callback_early_stopping(patience = 5)

  # Fit the model and save in history
  history <- # autoencoder %>%
    keras::fit(autoencoder, X, Y, validation_data = list(X, Y), loss = "mse",
               epochs = epochs, batch_size = batch_size, callbacks = list(checkpoint, early_stopping),
               verbose = ifelse(verbose,2,0), view_metrics = ifelse(verbose,"auto",0))


  # Final embeding
  embeddingSpace = stats::predict(encoder, Xtest)
  return(embeddingSpace)}
}

#' Preparing the input and output of the EDNN for a multiplex graph
#'
#' @param edge.list edge list as a dataframe with two columns.
#' @param edge.weight edge weights as a dataframe. Each column corresponds to a graph. By default, the \code{colnames} are considered as outcomes unless indicated in \code{outcome} argument.
#' @param outcome clinical outcomes for each graph. If not mentioned, the \code{colnames(edge.weight)} are considered by default.
#' @param indv.index the index of individual networks.
#' @param edge.threshold numeric value to set edge weights below the threshold to zero (default: 0). the greater edge weights do not change.
#' @param walk.rep number of repeats for the random walk (default: 100).
#' @param n.steps number of the random walk steps (default: 5).
#' @param random.walk boolean value to enable the random walk algorithm (default: TRUE).
#' @param verbose if \emph{TRUE} a progress bar is shown.
#'
#' @return the input and output required to train the EDNN
#' @export
#'
#' @examples
#' myNet = network_gen(N_nodes = 50)
#' graphData = myNet[["data_graph"]]
#' edge.list = graphData[,1:2]
#' edge.weight = graphData[,3:4]
#' XY = ednn_IOprepare(edge.list, edge.weight)
#' X = XY[["X"]]
#' Y = XY[["Y"]]
#'
ednn_IOprepare = function(edge.list, edge.weight, outcome=NULL, indv.index = NULL,
                          edge.threshold=0, walk.rep=10, n.steps=5,
                          random.walk = TRUE, verbose = TRUE){

  if (is.null(outcome))
    outcome = colnames(edge.weight)

  edge.list = t(edge.list)
  graph = igraph::graph(edge.list, directed = FALSE)

  N_nodes = length(igraph::V(graph))
  N_graph = ncol(edge.weight)

  if (is.null(indv.index))
    indv.index = 1:N_graph

  X = c()
  Y = c()
  outcome_node = c()
  individual_node = c()

  for (i in 1:N_graph){
    ## Set layer-specific weights for each graph and
    ## [optional] perform a thresholding if needed
    W_i = as.numeric(edge.weight[,i])
    graph_i = igraph::simplify(igraph::set.edge.attribute(graph, "weight", index=igraph::E(graph), W_i))
    graph_i = igraph::delete_edges(graph_i, igraph::E(graph_i)[igraph::E(graph_i)$weight < edge.threshold])

    ## Step 1) Input: Adjacency matrix calculation
    Adj_i = as.matrix(igraph::as_adj(graph_i,  attr = "weight"))

    ## Step 2) Output: Perform the fixed-length random walk
    ## and calculating the node visit probabilities.
    # Two options exist:
    # 1.repetitive simple random walks
    # 2.repetitive weighted random walks (specific to weighted graphs)
    if (random.walk){
      RW = RepRandomWalk (graph_i, Nrep = walk.rep, Nstep = n.steps, weighted_walk = TRUE, verbose = verbose)

      ## Step 3) Make it multilayer for EDNN
      X = rbind(X, Adj_i)
      Y = rbind(Y, RW$Probabilities)
    }else{
      X = rbind(X, Adj_i)
      Y = rbind(Y, Adj_i)
    }

    outcome_node = c(outcome_node, rep(outcome[i], N_nodes))
    individual_node = c(individual_node, rep(indv.index[i], N_nodes))
  }

  X = X / (apply(X, 1, sum) + .000000001)

  Result = list()
  Result[["X"]] = X
  Result[["Y"]] = Y
  Result[["outcome_node"]] = outcome_node
  Result[["individual_node"]] = individual_node
  return(Result)
}
