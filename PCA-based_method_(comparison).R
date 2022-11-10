
## Compare With Laplacian
graph = graph(t(NodeList), directed = FALSE)

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
  dist3[i] = (Q1[i,2:nPC] %*% Q2[i,2:nPC]) / 
  (sqrt(sum(Q1[i,2:nPC]^2))*sqrt(sum(Q2[i,2:nPC]^2)) + .000000001)
#dist3[i] = sum((Q1[i,2:nPC] - Q2[i,2:nPC])^2)

varNodes3 = order(dist3)[1:length(varNodes1)]

jc1 = length(intersect(varNodes,varNodes1))/length(union(varNodes,varNodes1))
jc2 = length(intersect(varNodes2,varNodes1))/length(union(varNodes2,varNodes1))
jc3 = length(intersect(varNodes3,varNodes1))/length(union(varNodes3,varNodes1))

JI1 = c(JI1, jc1)
JI2 = c(JI2, jc2)
JI3 = c(JI3, jc3)


l = length(NoiseSet)
NoiseSet[1] = 0
plot(NoiseSet[1:(l-1)], JI1[1:(l-1)], pch = 20, col = "red", type="b", 
     xlab="noise level", ylab="Jaccard index", log="x",
     ylim = c(0,1))
lines(NoiseSet[1:(l-1)], JI3[(l-1):1], pch=18, col="blue", type="b", xlab="x", ylab="y", log="x")
