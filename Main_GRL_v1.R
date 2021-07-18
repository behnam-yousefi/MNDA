## Graph Representation Learning
## By: Behnam Yousefi

rm(list = ls())
setwd("~/Desktop/R_Root/Microbiome/")
library(igraph)
library(keras)

RWR = function(G, start, MaxItt = 10, MaxSteps = 5, gamma = .5){
  Steps = MaxSteps
  RW_steps = c()
  for (i in 1:MaxItt)
    RW_steps = c(RW_steps, random_walk(G, start = node, steps = Steps))
  RW_steps = unique(RW_steps)
  return(RW_steps)
}

IterativeRandomWalk = function(G, Nrep = 10, Nstep = 5){
  N = length(V(G))
  X = matrix(0, (Nrep*N), N)
  Y = matrix(0, (Nrep*N), N)
  Xtest = c()
  
  i = 0
  for (node in V(G)){
    Neighbors = c(node,neighbors(G,node))
    NeighborsOneHot = rep(0,N)
    NeighborsOneHot[Neighbors] = 1
    Xtest = rbind(Xtest, NeighborsOneHot)
    
    for (rep in 1:Nrep){
      
      i = i + 1
      X[i,] = NeighborsOneHot
      # RW_steps = random_walk(G, start = node, steps = Nstep)
      RW_steps = RWR(G, start = node, MaxItt = 10, MaxSteps = Nstep)
      
      WalkOneHot = rep(0,length(V(G)))
      WalkOneHot[RW_steps] = 1
      Y[i,] = WalkOneHot
    }
  }
  
  colnames(X) = V(G)$names
  colnames(Y) = V(G)$names
  colnames(Xtest) = V(G)$names
  
  Result = list()
  Result[["X"]] = X
  Result[["Y"]] = Y
  Result[["Xtest"]] = Xtest
  
  return(Result)
}

## Dummy Data
G1 = sample_k_regular(10,5)
V(G1)$names = paste0("V",as.character(1:10))
plot(G1)
G2 = sample_k_regular(10,5)
V(G2)$names = paste0("V",as.character(11:20))
plot(G2)

## Read Microbiome Data
data = data.frame(read.delim(file="Data/Resulting_net_from_corr_eliminate_6M.txt",sep = " " ))
Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
G1 = graph(Nodelist, directed = FALSE)
IndvGraphWeights = data[,-1]
IndvGraphWeights = log(abs(IndvGraphWeights)+1)
IndvGraphWeights = (IndvGraphWeights - min(IndvGraphWeights)) / (max(IndvGraphWeights) - min(IndvGraphWeights))
AggGraphWeights = apply(IndvGraphWeights, 1, mean)
G1 = simplify(set.edge.attribute(G1, "weight", index=E(G1), AggGraphWeights))
G1 = delete_edges(G1, E(G1)[E(G1)$weight < .3])

data = data.frame(read.delim(file="Data/Resulting_net_from_corr_eliminate_9M.txt",sep = " " ))
Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
G2 = graph(Nodelist, directed = FALSE)
IndvGraphWeights = data[,-1]
IndvGraphWeights = log(abs(IndvGraphWeights)+1)
IndvGraphWeights = (IndvGraphWeights - min(IndvGraphWeights)) / (max(IndvGraphWeights) - min(IndvGraphWeights))
AggGraphWeights = apply(IndvGraphWeights, 1, mean)
G2 = simplify(set.edge.attribute(G2, "weight", index=E(G2), AggGraphWeights))
G2 = delete_edges(G2, E(G2)[E(G2)$weight < .3])

## Walker
RW1 = IterativeRandomWalk (G1, Nrep = 1, Nstep = 20)
RW2 = IterativeRandomWalk (G2, Nrep = 1, Nstep = 20)

# # Method I
# X = rbind(
#   cbind(RW1$X, rep(1,nrow(RW1$X)), rep(0,nrow(RW1$X))),
#   cbind(RW2$X, rep(0,nrow(RW2$X)), rep(1,nrow(RW2$X))) )
# Y = rbind(RW1$Y, RW2$Y)
# Xtest = rbind(
#   cbind(RW1$Xtest, rep(1,nrow(RW1$Xtest)), rep(0,nrow(RW1$Xtest))),
#   cbind(RW2$Xtest, rep(0,nrow(RW2$Xtest)), rep(1,nrow(RW2$Xtest))) )
# 
# # Method II
# X = rbind(
#   cbind(RW1$X, matrix(0,nrow(RW2$X), ncol(RW2$X))),
#   cbind(matrix(0,nrow(RW1$X), ncol(RW1$X)), RW2$X))
# Y = rbind(RW1$Y, RW2$Y)
# Xtest = rbind(
#   cbind(RW1$Xtest, matrix(0,nrow(RW2$Xtest), ncol(RW2$Xtest))),
#   cbind(matrix(0,nrow(RW1$Xtest), ncol(RW1$Xtest)), RW2$Xtest))

# Method III
X = rbind(RW1$X, RW2$X)
Y = rbind(RW1$Y, RW2$Y)
Xtest = rbind(RW1$Xtest, RW2$Xtest)


colnames(X) = paste0("V",as.character(1:ncol(X)))
colnames(Y) = paste0("V",as.character(1:ncol(Y)))
colnames(Xtest) = paste0("V",as.character(1:ncol(X)))


Nnode = ncol(X)
## Autoencoder traning

latentSize = 2
inputSize = ncol(X)
outputSize = ncol(Y)

{
# Define Encoder
enc_input = layer_input(shape = inputSize)
enc_output = enc_input %>% 
    layer_dense(units=latentSize, activation = "relu", 
                #kernel_regularizer = regularizer_l1(0.0001),
                #kernel_initializer = "Orthogonal", bias_initializer = "Zeros"
                #kernel_initializer = "Zeros", bias_initializer = "Zeros"
                )
#activation_tanh()
  
encoder = keras_model(enc_input, enc_output)
  
# Define decoder
dec_input = layer_input(shape = latentSize)
dec_output = dec_input %>% 
    layer_dense(units = outputSize, activation = "softmax") 
# layer_activation_relu()
  
decoder = keras_model(dec_input, dec_output)
# Define Auto-Encoder
aen_input = layer_input(shape = inputSize)
aen_output = aen_input %>% 
    encoder() %>% 
    decoder()
  
autoencoder = keras_model(aen_input, aen_output)
  
# Training configuration
autoencoder %>% compile(loss = "binary_crossentropy", optimizer = 'adam')
checkpoint <- callback_model_checkpoint(filepath = "My_model_temp", save_best_only = TRUE, verbose = 0)
early_stopping <- callback_early_stopping(patience = 5)
  
# Fit the model and save in history
#X = as.matrix(Prob_nodeVisit)
history <- autoencoder %>% fit(X, Y, validation_data = list(X, Y),
                                 epochs = 10, batch_size = 10, callbacks = list(checkpoint, early_stopping))

}

# Final embeding
latentSpace = encoder %>% predict(Xtest)
N = nrow(latentSpace) / 2
plot(latentSpace[,c(1,2)], cex = .2, pch = 20, col = c(rep(1,N),rep(2,N)))
# text(latentSpace[,c(1,2)],labels = 1:1144, col="red", cex = .7)

latentSpace_1 = latentSpace[1:N,]
latentSpace_2 = latentSpace[(N+1):(2*N),]

dist = rep(0,N)
G1 = set.vertex.attribute(G1, name = "name", value = names(V(G2)))
names(dist) = names(V(G1))
for (i in 1:N){
  #dist[i] = sum((latentSpace_1[i,] - latentSpace_2[i,])^2)
  dist[i] = (latentSpace_1[i,] %*% latentSpace_2[i,]) / ((sum(latentSpace_1[i,]^2))*(sum(latentSpace_2[i,]^2)))
}
hist(dist,50)

JI = rep(0,N)
names(JI) = names(V(G1))
for (i in 1:N){
  a = as.character(neighbors(G1, names(dist)[i]))
  b = as.character(neighbors(G2, names(dist)[i]))
  JI[i] = length(intersect(a,b))/length(union(a,b))
}
hist(JI,50)

cor(dist,JI)
plot(dist,JI,pch=20,cex=.5)

JI2 = rep(0,N)
names(JI2) = names(V(G1))
for (i in 1:N){
  n = neighbors(G1, names(dist)[i])
  for (j in n)
    n = c(n, neighbors(G1, j))
  a = as.character(unique(n))
  
  n = neighbors(G2, names(dist)[i])
  for (j in n)
    n = c(n, neighbors(G2, j))
  b = as.character(unique(n))
  
  JI2[i] = length(intersect(a,b))/length(union(a,b))
}
hist(JI2,50)

cor(dist,JI2)
plot(dist,JI2,pch=20,cex=.5)
