# Multiplex Network Differential Analysis (MNDA) <img width="100" alt="image" src="https://github.com/behnam-yousefi/MNDA/blob/master/Figures/MNDA_logo.png?raw=true">

Interactions between different biological elements are crucial for the function of biological systems and are usually represented as networks. The dysregulation of these networks can be associated with different clinical conditions such as diseases and response to treatments. In this work,  we propose an R package, called *multiplex network differential analysis* (MNDA) to quantify and test the variations in the local neighborhood of nodes between the two given conditions (e.g. case-control). We further show examples of finding important subnetworks in gene co-expression networks for response to treatment as well as an use case for individual specific networks (ISNs).

The core of the MNDA tool consists of three steps (Figure 1):
1. representing the nodes of all networks layers into a common embedding space (using EDNN);
2. calculate the distance between the nodes corresponding to the same element (e.g. gene);
3. detect the nodes whose neighborhood varies significantly based on statistical testing (using permuted graphs).

The current MNDA pipeline is designed for two conditions:\
a. two-layer network case corresponding to two (paired/unpaired) conditions (e.g. healthy-disease);\
b. multi-layer network case (e.g. ISNs) with two matched groups (e.g. before treatment-after treatment).

The future features to be considered:\
c. multi-layer network case (e.g. ISNs) with two unmatched groups (e.g. healthy-disease);\
d. multi-layer network case each layer with different condition (e.g. temporal networks).


**Note:** the first step (EDNN) is the same for all the cases but they are different in step 2 and 3.

**Note:** the network structures for different layers should be the same but different in weights.

Reference: [to be added]

## Installation
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

## Apply on simulated networks
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
the ```mnda_embedding_2layer()``` function represents all the nodes in a common embedding space (step 1); and the ```mnda_node_detection_2layer()``` function calculates the node-pair distances and assines a p-value to each node-pair (step 2 and 3). This process is repeated ```train.rep``` times to improve the robustness.

## Usage Example 1: drug response  
*MNDA pipeline for condition "a"*

In this example, we construct gene coexpression networks (GCNs) for drug responders and non-responders. We load gene expression profiles of cell lines of lung cancer as ```X``` and a binary vector of their response to the Tamoxifen drug as ```y```.
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

## Usage Example 2: application on individual specific networks
*MNDA pipeline for condition "b"*

In this analysis we consider a set of paired ISNs for two conditions (e.g. before treatment-after treatment) and a set of external variables for each individual (e.g. drug response and sex). The aim is to find nodes whose neighborhood variation between the two conditions is associated with the external variables.

In our example, we use the gene expression profile of blood samples before and after being stimulated with BCG vaccine and E coli. We then find genes whose neighborhood changes (dynamics) have a significant association with their sex. To impute ISNs, we use *lionessR* R package. We first read the ISN data to create the node list.
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

## Outreach
Any suggestions, collaboration and bug reports are welcome and appreciated.
Contact me via

Email: yousefi.bme@gmail.com,


[LinedIn](https://www.linkedin.com/in/behnam-yousefi-bme) or
[Twitter](https://twitter.com/behnam_bme)


