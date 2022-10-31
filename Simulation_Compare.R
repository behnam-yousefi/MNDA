## Simulation
## Graph Representation Learning
## By: Behnam Yousefi

rm(list = ls())
setwd("~/Desktop/R_Root/Microbiome_final/")
source('~/Desktop/R_Root/Microbiome_final/MNDA_functions.R')

library(igraph)
library(MASS)

## Read Microbiome Data
data = data.frame(read.delim(file="Data/Resulting_net_notNULL_MAGMACONF6M.txt",sep = " " ))
Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
MianGraph1 = graph(Nodelist, directed = FALSE)
IndvGraphWeights1 = data[,-1]
IndvGraphWeights1 = log(abs(IndvGraphWeights1)+1)
IndvGraphWeights1 = (IndvGraphWeights1 - min(IndvGraphWeights1)) / (max(IndvGraphWeights1) - min(IndvGraphWeights1))
AggGraphWeights1 = apply(IndvGraphWeights1, 1, mean)

rm(data, Nodelist)
N = length(V(MianGraph1))

JI1 = c()
JI2 = c()
JI3 = c()

NoiseSet = c(0, exp(-10:-2))
for (noise in NoiseSet){

  Threshold = 0
  G1 = simplify(set.edge.attribute(MianGraph1, "weight", index=E(MianGraph1), AggGraphWeights1))
  G1 = delete_edges(G1, E(G1)[E(G1)$weight < Threshold])
  A1 = as.matrix(as_adj(G1,  attr = "weight"))
  
  A2 = A1
  A2[1,] = A2[1,95:1]
  A2[2,] = A2[2,95:1]
  
  Aa1 = A1
  Aa2 = A2
  Gg2 = simplify(graph_from_adjacency_matrix(Aa2,  mode = "undirected", weighted = TRUE))
  Aa2 = as.matrix(as_adj(Gg2,  attr = "weight"))
  
  A1 = A1 + matrix(runif(9025 ,min = 0, max = noise), 95,95)
  A2 = A2 + matrix(runif(9025 ,min = 0, max = noise), 95,95)
  
  G1 = simplify(graph_from_adjacency_matrix(A1,  mode = "undirected", weighted = TRUE))
  G2 = simplify(graph_from_adjacency_matrix(A2,  mode = "undirected", weighted = TRUE))
  
  ## Walker
  RW1 = RandomWalkRestart (G1, Nrep = 100, Nstep = 5, weighted_walk = TRUE)
  RW2 = RandomWalkRestart (G2, Nrep = 100, Nstep = 5, weighted_walk = TRUE)

  A1 = as.matrix(as_adj(G1,  attr = "weight"))
  A2 = as.matrix(as_adj(G2,  attr = "weight"))

  X = rbind(A1, A2)
  X = X / (apply(X, 1, sum) + .000000001)
  Y = rbind(RW1$P, RW2$P)
  Xtest = X

  colnames(X) = paste0("V",as.character(1:ncol(X)))
  colnames(Y) = paste0("V",as.character(1:ncol(Y)))
  colnames(Xtest) = paste0("V",as.character(1:ncol(X)))


  Nnode = ncol(X)

  latentSize = 2
  inputSize = ncol(X)
  outputSize = ncol(Y)

  {
  # Define Encoder
  enc_input = layer_input(shape = inputSize)
  enc_output = enc_input %>% 
    layer_dense(units=latentSize, activation = "relu", 
                kernel_regularizer = regularizer_l2(0.0001),
                #kernel_initializer = "Orthogonal", bias_initializer = "Zeros"
                #kernel_initializer = "Zeros", bias_initializer = "Zeros"
    )
  #activation_tanh()
  
  encoder = keras_model(enc_input, enc_output)
  
  # Define decoder
  dec_input = layer_input(shape = latentSize)
  dec_output = dec_input %>% 
    layer_dense(units = outputSize, activation = "sigmoid")   # softmax
  
  decoder = keras_model(dec_input, dec_output)
  # Define Auto-Encoder
  aen_input = layer_input(shape = inputSize)
  aen_output = aen_input %>% 
    encoder() %>% 
    decoder()
  
  autoencoder = keras_model(aen_input, aen_output)
  
  # Training configuration binary_crossentropy
  autoencoder %>% compile(loss = "mse", optimizer = 'adam')
  checkpoint <- callback_model_checkpoint(filepath = "My_model_temp", save_best_only = TRUE, verbose = 0)
  early_stopping <- callback_early_stopping(patience = 5)
  
  # Fit the model and save in h
  history <- autoencoder %>% fit(X, Y, validation_data = list(X, Y), loss = "mse",
                                 epochs = 50, batch_size = 20, callbacks = list(checkpoint, early_stopping))
  
}

  # Final embeding
  latentSpace = encoder %>% predict(Xtest)
  
  latentSpace_1 = latentSpace[1:N,]
  latentSpace_2 = latentSpace[(N+1):(2*N),]
  Max = apply(latentSpace, 2, max)
  Min = apply(latentSpace, 2, min)
  
  plot(latentSpace_1, cex = .2, pch = 20, col = "red", xlim = c(Min[1],Max[1]), ylim = c(Min[2],Max[2]))
  text(latentSpace_1, labels = V(G1), cex = .6, col = "red", xlim = c(Min[2],Max[1]), ylim = c(Min[2],Max[2]))
  par(new=TRUE)
  plot(latentSpace_2, cex = .2, pch = 20, col= "blue", xlim = c(Min[1],Max[1]), ylim = c(Min[2],Max[2]))
  text(latentSpace_2, labels = V(G1), cex = .5, col= "blue", xlim = c(Min[2],Max[1]), ylim = c(Min[2],Max[2]))
  
  ## set with ground truth
  dist1 = rep(0,N)
  G1 = set.vertex.attribute(G1, name = "name", value = names(V(G2)))
  names(dist1) = names(V(G1))
  for (i in 1:N)
    dist1[i] = sum((Aa1[i,] - Aa2[i,])^2)
  
  varNodes1 = which(dist1>0)
  
  ## Compare with our method
  dist = rep(0,N)
  G1 = set.vertex.attribute(G1, name = "name", value = names(V(G2)))
  names(dist) = names(V(G1))
  for (i in 1:N)
    dist[i] = (latentSpace_1[i,] %*% latentSpace_2[i,]) / (sqrt(sum(latentSpace_1[i,]^2))*sqrt(sum(latentSpace_2[i,]^2)) + .000000001)
  # dist[i] = sum((latentSpace_1[i,] - latentSpace_2[i,])^2)
  # dist[i] = (latentSpace_1[i,] %*% latentSpace_2[i,])
  varNodes = order(dist)[1:length(varNodes1)]
  
  ## Compare With Laplacian
  L1 = as.matrix(laplacian_matrix(G1, normalized = TRUE))
  L2 = as.matrix(laplacian_matrix(G2, normalized = TRUE))
  L1[is.na(L1)] = 0
  L2[is.na(L2)] = 0
  
  dist2 = rep(0,N)
  G1 = set.vertex.attribute(G1, name = "name", value = names(V(G2)))
  names(dist2) = names(V(G1))
  for (i in 1:N)
    dist2[i] = (L1[i,] %*% L2[i,]) / (sqrt(sum(L1[i,]^2))*sqrt(sum(L2[i,]^2)) + .000000001)
  
  varNodes2 = order(dist2)[1:length(varNodes1)]
  
  ## Compare With PCA of Laplacian
  PCA = prcomp(L1, scale. = FALSE)
  Q1 = PCA$rotation[,ncol(PCA$rotation):1]
  PCA = prcomp(L2, scale. = FALSE)
  Q2 = PCA$rotation[,ncol(PCA$rotation):1]
  
  nPC = latentSize + 1
  dist3 = rep(0,N)
  G1 = set.vertex.attribute(G1, name = "name", value = names(V(G2)))
  names(dist3) = names(V(G1))
  for (i in 1:N)
    dist3[i] = (Q1[i,2:nPC] %*% Q2[i,2:nPC]) / (sqrt(sum(Q1[i,2:nPC]^2))*sqrt(sum(Q2[i,2:nPC]^2)) + .000000001)
  #dist3[i] = sum((Q1[i,2:nPC] - Q2[i,2:nPC])^2)
  
  varNodes3 = order(dist3)[1:length(varNodes1)]
  
  jc1 = length(intersect(varNodes,varNodes1))/length(union(varNodes,varNodes1))
  jc2 = length(intersect(varNodes2,varNodes1))/length(union(varNodes2,varNodes1))
  jc3 = length(intersect(varNodes3,varNodes1))/length(union(varNodes3,varNodes1))
    
  JI1 = c(JI1, jc1)
  JI2 = c(JI2, jc2)
  JI3 = c(JI3, jc3)
}

l = length(NoiseSet)
NoiseSet[1] = 0
plot(NoiseSet[1:(l-1)], JI1[1:(l-1)], pch = 20, col = "red", type="b", 
     xlab="noise level", ylab="Jaccard index", log="x",
     ylim = c(0,1))
lines(NoiseSet[1:(l-1)], JI3[(l-1):1], pch=18, col="blue", type="b", xlab="x", ylab="y", log="x")
