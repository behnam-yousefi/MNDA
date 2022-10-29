## Multiplex Graph Representation Learning and
## multiplex network differential analysis
## All ISNs in one embedding
## Graph Representation Learning
## By: Behnam Yousefi

library(igraph)
library(keras)

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

RepRandomWalk = function(graph, Nrep = 10, Nstep = 5, weighted_walk = TRUE){
  # Repetitive Fixed-length [weighted] random walk algorithm
  
  N = length(V(graph))
  X = matrix(0, N, N)
  Y = matrix(0, N, N)
  
  i = 0
  for (node in V(graph)){
    Neighbors = c(node,neighbors(graph,node))
    NeighborsOneHot = rep(0,N)
    NeighborsOneHot[Neighbors] = 1
    
    i = i+1
    X[i,] = NeighborsOneHot
    
    WalkRep = matrix(0, Nrep, N)
    for (rep in 1:Nrep){
      
      if (weighted_walk)
        RW_steps = WeightdRandomWalk(graph, startNode = node, maxStep = Nstep)
      else
        RW_steps = random_walk(graph, start = node, steps = Nstep)
      
      WalkRep[rep,RW_steps] = 1
    }
    
    Y[i,] = apply(WalkRep, 2, sum)
  }
  
  P = Y / Nrep
  P[P>1] = 1
  diag(P) = 1
  
  colnames(X) = V(graph)$names
  colnames(Y) = V(graph)$names
  colnames(P) = V(graph)$names
  
  Result = list()
  Result[["X"]] = X
  Result[["Y"]] = Y
  Result[["P"]] = P
  
  return(Result)
}

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

