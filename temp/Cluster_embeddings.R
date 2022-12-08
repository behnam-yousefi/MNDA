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

# The co-association similarity can be calculated as
CS = coCluster_matrix(Clusters)
pheatmap::pheatmap(CS)

# Then a PCC clustering approach can be used for clustering along with obtaining the 
# number of clusters
CM = consensus_matrix(CS, max.cluster = 5, resample.ratio = 0.7, 
                      max.itter = 50, clustering.method = "pam")
Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, pch = 20, type = "b", ylab = "robustness score", xlab = "number of clusters")

Kopt = Scores[["Kopt_LogitScore"]]
clusters = PamClustFromAdjMat (CM[[Kopt]], k = 2, alpha = 1, adj.conv = FALSE)


## Co-association similarity measure between the corresponding nodes (in a 2-layer network case)
N_nodes = ncol(CS)/2
Sim_crsp_nd = c()
for (i in 1:N_nodes)
  Sim_crsp_nd = c(Dist_crsp_nd, CS[N_nodes + i])

plot(sort(Sim_crsp_nd), pch = 20)
high_var_nodes = order(Sim_crsp_nd, decreasing = FALSE)[1:9]


