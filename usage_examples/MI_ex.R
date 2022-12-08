rm(list=ls())

setwd("~/Desktop/R_Root/MNDA/Applications/")
data = read.table("~/Desktop/R_Root/SLEmap/Data/MI/NANOSTRING_7COND.database.txt", header = TRUE)
Pheno = data[,c(1:7)]
y = Pheno$STIMULUS_NAME
X = as.matrix(data[,-c(1:7,568)])

adj_cont = abs(cor(X[y=="Null",]))
adj_case = abs(cor(X[y=="BCG",]))

adj_list = list(adj_cont, adj_case)
graph_data = as.mnda.graph(adj_list, outcome = c("cont","case"))

embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep=50, walk.rep=10,
                                           epochs = 10, batch.size = 20,
                                           random.walk=FALSE, null.perm = FALSE)

Results = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "BH")
Nodes = Results$significant_nodes
# write.table(Nodes, "~/Desktop/nodes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

NodeSet = Nodes

graph_to_plot = cbind(graph_data[,1:2], W = graph_data[,3] - graph_data[,4])
hist(graph_to_plot$W)

G = mnda::as.igraph(graph_to_plot, 1.0)

pdf("~/Desktop/BCG.pdf")
subgraph_plot(G, NodeSet)
dev.off()

