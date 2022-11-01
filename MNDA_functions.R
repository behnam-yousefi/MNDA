## Multiplex Graph Representation Learning and
## multiplex network differential analysis
## By: Behnam Yousefi

library(igraph)
library(keras)

### Weighted Random Walk algorithm
WeightdRandomWalk = function(graph, startNode, maxStep = 5, node_names = FALSE){
  # graph: An igraph object
  # startNode: the starting node can be name or index.
  # node_names ???
  
  # output value: The set of nodes passed by the random walker.
  
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


### Repetitive Fixed-length [weighted] random walk algorithm
RepRandomWalk = function(graph, Nrep = 10, Nstep = 5, weighted_walk = TRUE){
  # graph: An igraph object
  
  # output values:
  # Steps (S): The total number of times a node is visited starting from the corresponding node in the row.
  # Probabilities (P): The node visit probabilities starting from the corresponding node in the row.
  
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


### Encoder decoder neural network (EDNN) function
EDNN = function(X, Y, Xtest, latentSize = 2, epochs = 10, batch_size = 100, l2reg = 0){

  Nnode = ncol(X)
  inputSize = ncol(X)
  outputSize = ncol(Y)

  # Define Encoder
  enc_input = layer_input(shape = inputSize)
  enc_output = enc_input %>% 
    layer_dense(units=latentSize, activation = "relu",
                kernel_regularizer = regularizer_l2(l2reg))
  
  encoder = keras_model(enc_input, enc_output)
  
  # Define decoder
  dec_input = layer_input(shape = latentSize)
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
  latentSpace = encoder %>% predict(Xtest)
  return(latentSpace)
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

### function to calculate rank
Rank = function(x, decreasing = FALSE) {
  # hint: What is the difference between Order and Rank
  # Order: [the index of the greatest number, ..., the index of the smallest number]
  # Rank: [the rank of the 1st number, ..., the rank of the last number]
  # In Rank, the order of the numbers remains constant,so can be used for ranksum.
  # ex) > a = c(10, 20, 50, 30, 40)
  #     > order(a)
  #     [1] 1 2 4 5 3
  #     > Rank(a)
  #     [1] 1 2 5 3 4
  
  Order = order(x, decreasing = decreasing)
  
  Rank = rep(0,length(Order))
  for(i in 1:length(Order))
    Rank[Order[i]] = i
  
  return(Rank)
}

