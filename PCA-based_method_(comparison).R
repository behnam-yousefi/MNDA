
graph = graph(t(NodeList), directed = FALSE)

W_1 = as.numeric(EdgeWeights[,1])
W_2 = as.numeric(EdgeWeights[,2])

graph_1 = simplify(set.edge.attribute(graph, "weight", index=E(graph), W_1))
graph_2 = simplify(set.edge.attribute(graph, "weight", index=E(graph), W_2))

L1 = as.matrix(laplacian_matrix(graph_1, normalized = TRUE))
L2 = as.matrix(laplacian_matrix(graph_2, normalized = TRUE))
L1[is.na(L1)] = 0
L2[is.na(L2)] = 0

## Laplacian distance
Dist_LapMat = rep(0,N_nodes)
for (i in 1:N_nodes)
  Dist_LapMat[i] = (L1[i,] %*% L2[i,]) / (sqrt(sum(L1[i,]^2))*sqrt(sum(L2[i,]^2)) + .000000001)

plot(sort(Dist_LapMat))
order(Dist_LapMat, decreasing = FALSE)[1:5]

## PCA-based distance
PCA = prcomp(L1, scale. = FALSE)
Q1 = PCA$rotation[,ncol(PCA$rotation):1]
PCA = prcomp(L2, scale. = FALSE)
Q2 = PCA$rotation[,ncol(PCA$rotation):1]

nPC = 6
Dist_PCA = rep(0,N_nodes)
for (i in 1:N_nodes)
  Dist_PCA[i] = (Q1[i,2:nPC] %*% Q2[i,2:nPC]) / (sqrt(sum(Q1[i,2:nPC]^2))*sqrt(sum(Q2[i,2:nPC]^2)) + .000000001)
  # Dist_PCA[i] = sum((Q1[i,2:nPC] - Q2[i,2:nPC])^2)

plot(sort(Dist_PCA))
order(Dist_PCA, decreasing = FALSE)[1:5]


# G1 = set.vertex.attribute(G1, name = "name", value = names(V(G2)))
# names(dist2) = names(V(G1))
