# Supplementary File

## 1. Methods overview

MNDA is a computational tool for Multiplex Network Differential Analysis that operates on multiplex networks to detect nodes whose neighborhoods have significant variations across layers. The core of the tool consists of three steps:

1.	representing the nodes of all network layers into a common embedding space (using an encoder-decoder neural network - EDNN);
2.	calculating the distance between multiplex corresponding nodes (e.g. genes in biological networks);
3.	detecting the nodes whose neighborhoods vary significantly based on statistical testing (using permuted graphs).

The EDNN is composed of shallow encoder-decoder neural networks with the number of inputs and outputs being equal to the number of nodes in one layer (the nodes are the same from one layer to the other). For each node, the encoder input is a vector of its connection weights with the other nodes. If no link exists between two nodes, the corresponding input is set to zero. The decoder output for each node is a vector of node visit probabilities calculated based on a repeated fixed-length weighted random walk algorithm (Yousefi et al., 2023). The bottleneck layer serves as an embedding space in which dissimilarity between corresponding nodes are computed. The default dissimilarity measure implemented in the *mnda* R package is the cosine-distance. For technical details and motivations for parameter settings related to steps 1 and 2, we refer to (Yousefi et al., 2023). 

Step 3. involves assessing the significance of the calculated dissimilarities between corresponding nodes. Currently, the MNDA accommodates the following biological network analysis scenarios:

a.	Given two groups of samples corresponding to two conditions (e.g. independent healthy-diseased, male-female) and group-level biological networks $N_{C1}(V_{C1}, E_{C1})$ and $N_{C2}(V_{C2}, E_{C2})$ for each condition, form pairs of corresponding nodes ( $n_{1j}$ , $n_{2j}$ ); $n_{ij}$ is a node of $N_{Ci}$ , $i=1,2$ and $j$ runs from $1$ to the cardinality of $V_{Ci}$. 

b.	Given a single group of samples, for which two states are available $S_1$ and $S_2$ (e.g. matched case-control pairs, pre-post treatment for the same sample) and sample-level biological networks $N_{k,S_1}$ and $N_{k,S_2}$ for each state and individual $k$, form pairs of corresponding nodes ( $n_{k,1j}$ , $n_{k,2j}$ ); $n_{k,ij}$ is a node of $N_{k,S_i}$ , $i=1,2$ and $j$ runs from $1$ to the cardinality of each network (all networks are assumed to have the same node cardinality). 

Whereas in context *a* only one vector of dissimilarities is created with a length corresponding to the number of nodes in the group-level biological networks, in context *b* – the scenario of sample-specific or individual-specific networks - multiple such vectors are created, one for each sample. In both contexts, an empirical null distribution of the chosen dissimilarity statistic is computed by reshuffling edges in the available networks (2 in scenario *a* and 2 times the number of matched sample pairs in scenario *b*). The user can specify whether to keep the node degree distribution intact for each network. Raw P-values thus obtained for each node annotation (corresponding nodes in a node pair have the same annotation) can or cannot be adjusted for multiple testing. For scenario *a* as many tests as the cardinality of $V_{C1}$  ( $V_{C2}$ ) are carried out. In context *b* the number of tests corresponds to the total number of nodes across samples for a single state. *Bonferroni* correction is the default option and the user can opt for less conservative correction methods.

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
This is required only once for the installation.

### 2.2. Apply on simulated networks

To test the ```mnda``` package, a toy example multilayer network can be generated using the ```network_gen()``` function:
`````{R}
myNet = network_gen(N_nodes = 100, N_var_nodes = 5, N_var_nei = 90, noise_sd = .01)
`````
The process is described as the following:
1. two fully connected networks with ```N_nodes``` number of matching nodes and uniform random edge weights is generated.
2. a number of ```N_var_nodes``` nodes are randomly selected to have different edge weights with ```N_var_nei``` number of nodes (```N_var_nei``` $<$ ```N_nodes``` $- 1$) between the two networks.
3. a set of random *Gaussian* noise with zero mean and sd = ```noise_sd``` is generated and added to all of the edge weights.

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
the ```mnda_embedding_2layer()``` function represents all the nodes in a common embedding space (step 1); and the ```mnda_node_detection_2layer()``` function calculates the node-pair distances and asignes a P-value to each node-pair (step 2 and 3). This process is repeated ```train.rep``` times to improve the robustness. The source code available at [usage_examples/network_generation_ex.R](https://github.com/behnam-yousefi/MNDA/blob/master/usage_examples/network_generation_ex.R).

### 2.3. Usage Example 1: drug response  

In this example, which is a showcase for condition *a*, we construct gene coexpression networks (GCNs) for drug responders and non-responders. To this end, we use the PRISM dataset (Corsello et al., 2020), which is a cell line-based drug screening dataset. To reduce the dimensionality, 2000 genes that are highly variant across all the cell lines are selected and reposited. The gene expression profile of lung cancer cell lines as ```X``` and a binary vector of their response to the *Tamoxifen* drug as ```y``` can be loaded accordingly:
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
* the random walk algorithm can be disabled by ```random.walk = FALSE``` to decrease the running time;
* the network permutation and representation can be disabled by ```null.perm = FALSE``` to decrease the running time;
* the calculated P-values can be adjusted by setting a method in the ```p.adjust.method``` argument.

### 2.4. Usage Example 2: application on individual specific networks

In this example we use the data of Milieu Interieur project (Thomas et al., 2015; Piasecka et al., 2018), where immune transcriptional profiles of bacterial-, fungal-, and viral- induced blood samples in an age- and sex- balanced cohort of 1,000 healthy individuals are generated. Here, the aim would be to find genes whose neighborhood significantly varies between the two conditions of stimulated and unstimulated. Following the MNDA pipeline, we first construct a set of paired ISNs for the two conditions, i.e. before and after stimulation, using the *lionessR* R package (Kuijjer et al., 2019 a; Kuijjer et al., 2019 b). In each network, nodes and edge weights represent genes and the correlation of their expressions, respectively. The imputed ISNs are reposited in ```"usage_examples/Data/ISN_net.rds"```. We first read the ISN data and create the node list.
`````{R}
data = data.frame(readRDS("Data/ISN_BCG.rds"))
nodeList = t(sapply(rownames(data), function(x) strsplit(x,"_")[[1]]))
`````
Next, we create the individual variable *data.frame* with three columns of *Individual IDs* (indecis), *stimulation condition*, and *sex* (F,M).
`````{R}
y = colnames(data)
y = data.frame(t(data.frame(strsplit(y, "_"))))
colnames(y) = c("ID", "Stim", "Sex")
`````
Now that we have all the ISNs with their phenotypes, we can perform two types of analysis on the population level:

1- Aggregate ISNs (by averaging) into two groups (i.e. pre/post stimulation) and find genes with significant neighbourhood variation.
This will be similar to *Usage Example 1* in context *a* (see above). We first obtain the two aggregated networks of pre- and post- stimulation;
`````{R}
data_agg = cbind(apply(data[,y$Stim=="Null"], 1, mean),
                 apply(data[,y$Stim=="BCG"], 1, mean))
graph_data = cbind(nodeList, data_agg)
colnames(graph_data) = c("V1", "V2", "Null", "Stim")
`````
and then perform the two-layer MNDA pipeline.
`````{R}
embeddingSpaceList = mnda_embedding_2layer(graph_data, edge.threshold = .1,
                                           train.rep = 50, epochs = 25, batch.size = 10,
                                           random.walk = FALSE, null.perm = FALSE)
mnda_output = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "bonferroni", alpha = .01)
`````

2- Project nodes of all the ISNs in the same embedding space and find significant genes in context *b* (see above).
In this analysis, the ISNs of pre- and post- stimulation should be paired. Therefore, for each individual-gene, we have two points in the embedding space: one correspond to pre-stimulation and the other correspond to post-stimulation. Calculating the distance between these pairs, we will have a matrix of distances of size $N_{individual} \times N_{gene}$.

To implement this, we use ```mnda_embedding()``` and ```mnda_node_distance()``` commands, respectively.
`````{R}
graph_data = cbind(nodeList, data)
embeddingSpaceList = mnda_embedding(graph_data, outcome = y$Stim, indv.index = y$ID,
                                    train.rep=50, walk.rep=10, epochs=10, batch.size=50,
                                    random.walk=FALSE)
Dist = mnda_node_distance(embeddingSpaceList)
`````
Having the distance matrix, one can find extreme distances, i.e. highly variable/constant gene neighbourhoods, or find association of any variable with them, i.e. genes whose neighbourhoods are significantly associated with a variable such as sex:
`````{R}
sex = y[duplicated(y$ID), "Sex"]
Pval = mnda_distance_test_isn(Dist, p.adjust.method = "bonferroni")
`````


The source code is available at "[usage_examples/](https://github.com/behnam-yousefi/MNDA/blob/master/usage_examples/)"

## References
Corsello,S.M. et al. (2020) Discovering the anticancer potential of non-oncology drugs by systematic viability profiling. Nature Cancer, 1, 235–248.\
Kuijjer,M.L. et al. (2019 a) Estimating Sample-Specific Regulatory Networks. iScience, 14, 226–240.\
Kuijjer,M.L. et al. (2019 b) lionessR: single sample network inference in R. BMC Cancer, 19, 1003.\
Piasecka,B. et al. (2018) Distinctive roles of age, sex, and genetics in shaping transcriptional variation of human immune responses to microbial challenges. Proc. Natl. Acad. Sci. U. S. A., 115, E488–E497.\
Thomas,S. et al. (2015) The Milieu Intérieur study—an integrative approach for study of human immunological variance. Clin. Immunol., 157, 277–293.\
Yousefi,B. et al. (2023) Capturing the dynamics of microbiomes using individual-specific networks. bioRxiv 2023.01.22.525058.

