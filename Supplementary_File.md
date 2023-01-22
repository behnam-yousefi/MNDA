# Supplementary File

## 1. Methods overview

Multiplex network differential analysis (MNDA) is a computational tool implemented for multiplex networks to detect nodes whose neighborhoods have significant variations. The core of the MNDA tool consists of three steps:

1. representing the nodes of all networks layers into a common embedding space (using EDNN);
2. calculate the distance between the nodes corresponding to the same element (e.g. gene);
3. detect the nodes whose neighborhood varies significantly based on statistical testing (using permuted graphs). 

The EDNN is composed of shallow encoder-decoder neural networks with the number of inputs and outputs being equal to the number of nodes in one layer (the nodes are the same from one layer to the other). For each node, the encoder input is a vector of its connection weights with the other nodes. If no link exists between two nodes, the corresponding input is set to zero. The decoder output for each node is a vector of node visit probabilities calculated based on a *repeated fixed-length weighted random walk algorithm* (Yousefi et al.). The bottleneck layer serves as an embedding space in which dissimilarity between corresponding nodes are computed. The default dissimilarity measure implemented in the MNDA R package is the cosine-distance. For technical details and motivations for parameter settings related to steps 1 and 2, we refer to (Yousefi et al.). Step 3. involves assessing the significance of the calculated dissimilarities between corresponding nodes. Currently, the MNDA package considers two  scenarios: 

a. comparison of two groups of independent samples (e.g. network for healthy and disease populations); \
b. comparison of multiple matched samples (e.g. multiple individual-specific networks, matched according to pre- and post-treatment administration). \
In the two-layer network case, each node pair corresponds to a distance measure and its significance is assessed on the basis of a null distribution. On the other hand, each node pair in the multi-layer network case corresponds to a set of distances across all the individuals that are classified into two groups. Therefore, for each node pair, the distances can be classified into two sets and a two sample test (e.g. t-test of Wilcoxon-test) can be used for significance assessment. 

## 2. Implementation in R
### 2.1. Installation
Install from CRAN
`````{R}
install.packages("mnda")
`````
Install the latest version from GitHub
`````{R}
devtools::install_github("behnam-yousefi/MNDA/package/mnda")
`````
MNDA will also install TensorFlow and Keras for R, which need to be activated by installation of Miniconda. For this, according to the installation guidline of Keras ([here](https://cran.r-project.org/web/packages/keras/vignettes/index.html)):
`````{R}
library(keras)
install_keras()
`````
This is only required once for the installation.

## 2.2. Apply on simulated networks
*MNDA pipeline for condition "a"*

To test the ```mnda``` package, a toy example multilayer network can be generated using the ```network_gen()``` function:
`````{R}
myNet = network_gen(N_nodes = 100, N_var_nodes = 5, N_var_nei = 90, noise_sd = .01)
`````
The process is described as the following:
1. two identical fully connected networks with ```N_nodes``` number of nodes and uniform random edge weights is generated.
2. a number of ```N_var_nodes``` nodes are randomly selected to have different edge weights with ```N_var_nei``` number of nodes (```N_var_nei``` $<$ ```N_nodes``` $- 1$) between the two networks.
3. a set of random Gaussian noise with zero mean and sd = ```noise_sd``` is generated and added to all of the edge weights.

The generated multiplex network and the set of the randomly selected nodes are accessible by the following lines, respectively.
`````{R}
graph_data = myNet$data_graph
var_nodes = myNet$var_nodes
`````
We then feed ```graph_data``` to the MNDA pipeline specialized for a two-layer multiplex network (condition *"a"*), which is composed of two commands:
`````{R}
embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep = 50)
mnda_output = mnda_node_detection_2layer(embeddingSpaceList)
print(mnda_output$high_var_nodes_index)
`````
the ```mnda_embedding_2layer()``` function represents all the nodes in a common embedding space (step 1); and the ```mnda_node_detection_2layer()``` function calculates the node-pair distances and assines a p-value to each node-pair (step 2 and 3). This process is repeated ```train.rep``` times to improve the robustness. The source code available at [usage_examples/network_generation_ex.R](https://github.com/behnam-yousefi/MNDA/blob/master/usage_examples/network_generation_ex.R).

## 2.3. Usage Example 1: drug response  
*MNDA pipeline for condition "a"*

In this example, we construct gene coexpression networks (GCNs) for drug responders and non-responders. To this end, we use the PRISM dataset (Corsello et al., 2020), which is a cell line-based drug screening dataset. To reduce the dimensionality, 2000 genes that are highly variant across all the cell lines are selected and reposited. The gene expression profile of lung cancer cell lines as ```X``` and a binary vector of their response to the Tamoxifen drug as ```y``` can be loaded accordingly:
`````{R}
data = readRDS("Data/GCN2Layer_data_lung_tamoxifen_2000genes.rds")
X = data[[1]]
y = data[[2]]
`````
Next, we construct adjacency matrices of GCN for each condition,
`````{R}
adj_res = abs(cor(X[y=="res",]))
adj_nonres = abs(cor(X[y=="non_res",]))
diag(adj_res) = 0
diag(adj_nonres) = 0
`````
and convert them to the ```mnda``` multiplex network format.
`````{R}
adj_list = list(adj_res, adj_nonres)
graph_data = as.mnda.graph(adj_list, outcome = c("res","non_res"))
`````
Now we can call the MNDA pipeline as in the previous example
`````{R}
embeddingSpaceList = mnda_embedding_2layer(graph_data, edge.threshold = .1, train.rep = 50, epochs = 20, batch.size = 10, random.walk = FALSE, null.perm = FALSE)
mnda_output = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "bonferroni")
Nodes = mnda_output$high_var_nodes
`````
**hints for large networks:** 
* the random walk algorithm can be disabled by ```random.walk = FALSE``` for the sake of running time;
* the network permutation and representation can be disabled by ```null.perm = FALSE``` for the sake of running time;
* the calculated p-values can be adjusted by setting a method in the ```p.adjust.method``` argument.

The source code available at [usage_examples/drug_response_ex.R](https://github.com/behnam-yousefi/MNDA/blob/master/usage_examples/drug_response_ex.R)

## 2.4. Usage Example 2: application on individual specific networks
*MNDA pipeline for condition "b"*

In this example we use the data of Milieu Interieur project (Thomas et al., 2015; Piasecka et al., 2018), where immune transcriptional profiles of bacterial-, fungal-, and viral-induced blood samples in an age- and sex-balanced cohort of 1,000 healthy individuals are generated. Here, the aim would be to find genes whose neighborhood variation, between the two conditions stimulated and unstimulated, is associated with sex. Following the MNDA pipeline, we first construct a set of paired ISNs for the two conditions, i.e before stimulation and after treatment using the *lionessR* R package (Marieke Lydia Kuijjer et al., 2019; Marieke L. Kuijjer et al., 2019). In each network, nodes represent genes and the edge weights demonstrate the correlation of gene expression. The imputed ISNs are reposited in ```"usage_examples/Data/ISN_net.rds"```. We first read the ISN data creat the node list.
`````{R}
data = data.frame(readRDS("Data/ISN_net.rds"))
nodeList = t(sapply(rownames(data), function(x) strsplit(x,"_")[[1]]))
`````
Next, we create the individual variable data frame with three columns: Individual IDs (indecis), stimulation condition, sex (F,M).
`````{R}
y = colnames(data)
y = data.frame(t(data.frame(strsplit(y, "_"))))
`````
We then form the ```graph_data``` and call ```mnda_embedding()``` function to embed genes into the embedding space.
`````{R}
graph_data = cbind(nodeList, data)
embeddingSpaceList = mnda_embedding(graph_data, outcome = y[,2], indv.index = y[,1],train.rep=2, walk.rep=10, epochs=5, batch.size=100,random.walk=FALSE)
`````
For each individual-gene pair, we calculate the embedding distance between the nodes corresponding to before and after stimulation using ```mnda_node_distance()```. This results in distance matrices of *individual-by-gene*. Finally, we test the association of sex with the distance of each gene using ```mnda_distance_test_isn()```.
`````{R}
Dist = mnda_node_distance(embeddingSpaceList)
Pval = mnda_distance_test_isn(Dist, rep(c(1,2),25))
`````

The source code available at [usage_examples/ISN_ex.R](https://github.com/behnam-yousefi/MNDA/blob/master/usage_examples/ISN_ex.R)


## References
Corsello,S.M. et al. (2020) Discovering the anticancer potential of non-oncology drugs by systematic viability profiling. Nature Cancer, 1, 235–248.\
Kuijjer,M.L. et al. (2019) Estimating Sample-Specific Regulatory Networks. iScience, 14, 226–240.\
Kuijjer,M.L. et al. (2019) lionessR: single sample network inference in R. BMC Cancer, 19, 1003.\
Piasecka,B. et al. (2018) Distinctive roles of age, sex, and genetics in shaping transcriptional variation of human immune responses to microbial challenges. Proc. Natl. Acad. Sci. U. S. A., 115, E488–E497.\
Thomas,S. et al. (2015) The Milieu Intérieur study—an integrative approach for study of human immunological variance. Clin. Immunol., 157, 277–293.


