## Obtaining the embedding space using EDNN (e.g. run MNDA_2LayerNet.R),
## where the nodes in all the layers of the network are mapped to a same space,
## one can be stratify these nodes into some clusters to analyse neighborhood transitions.

## The input to this script is a list of embedding spaces (e.g. embeddingSpaceList)
embeddingSpaceList = readRDS("Data/Embedding_Space/Embedding_Space_1.rds")

# remove the element "outcome" from embeddingSpaceList 
embeddingSpaceList_temp = embeddingSpaceList[-length(embeddingSpaceList)]

# Here we use a method based on method perturbation consensus clustering (MPCC) to
# [i] obtain a stable clustering given different embedding spaces
# [ii] obtain the optimum number of clusters
# [iii] (optional) obtain the co-association similarity measure which can also be used to
#       detect the highly and lowly variable nodes
source("~/Desktop/R_Root/ConsensusClustering/CC_functions.R")

# Considering each repeat of the embedding space as a data view,
# we can perform a multivitamin clustering:
Clusters = multiview_kmeans(embeddingSpaceList_temp, rep = 50, 
                            range.k = c(5,40), method = "random")
Adj = coCluster_matrix(Clusters)
pheatmap::pheatmap(Adj)

## Co-association similarity measure between the corresponding nodes (in a 2-layer network case)
N_nodes = ncol(Adj)/2
Sim_crsp_nd = c()
for (i in 1:N_nodes)
  Sim_crsp_nd = c(Dist_crsp_nd, Adj[N_nodes + i])

high_var_nodes = order(Sim_crsp_nd, decreasing = FALSE)[1:9]


