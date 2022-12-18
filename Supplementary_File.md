# Supplementary File

## 1 Methods overview

Multiplex network differential analysis (MNDA) is a computational tool implemented for multiplex networks to detect nodes whose neighborhoods have significant variations. The core of the MNDA tool consists of three steps:

1. representing the nodes of all networks layers into a common embedding space (using EDNN);
2. calculate the distance between the nodes corresponding to the same element (e.g. gene);
3. detect the nodes whose neighborhood varies significantly based on statistical testing (using permuted graphs). 

The EDNN is composed of shallow encoder-decoder neural networks with the number of inputs and outputs being equal to the number of nodes in one layer (the nodes are the same from one layer to the other). For each node, the encoder input is a vector of its connection weights with the other nodes. If no link exists between two nodes, the corresponding input is set to zero. The decoder output for each node is a vector of node visit probabilities calculated based on a *repeated fixed-length weighted random walk algorithm (see below). The bottleneck layer is finally used as the embedding space for all the nodes in both layers.

*Repeated fixed-length weighted random walk algorithm.* We define the node visit probabilities for the decoder output as the probability of a random walker passing node $j$ starting from node $i$, $P(i|j)$. This enables us to characterise the local structures of the networks. The implemented random walker has two properties: It is weighted and fixed-length. A weighted random walker considered the edge weights to choose each step. The probability of moving from node $i$ to node $j$ is proportional to the edge weight $w_{ij}$ linking node $i$ to node $j$:

$$ P(i\xrightarrow{}j)=\frac{w_{ij}}{\Sigma_{j\in{}\mathcal{N}_i}|w_{ij}|} $$

The walk length is also set to a constant value to keep it local around the node queried node. The random walk process is repeated for several times to calculate the node visit probabilities. Although the idea of weighted random walks is not new, to our knowledge, no customized code was available for their use in our framework. Our implementation of the random walk algorithm is available for the users in the package. 

After the EDNN is trained, the nodes are embedded into a low dimensional vector space, based on which we then calculate the distance between the corresponding node pairs, i.e nodes that correspond to the same object in two networks. By default, we use the *cosine* distance which is defined as

$$ d_{cos}(A,B)=1-\frac{A.B}{|A|.|B|} $$

For two vectors of $A$ and $B$. 



## 1.1. fixed-length weighted random walk algorithm

## 1.1. fixed-length weighted random walk algorithm

## 2. Implementation in R
### 2.1 Installation
Inatall from CRAN
`````{R}
install.packages("mnda")
`````
Inatall the latest version from GitHub
`````{R}
devtools::install_github("behnam-yousefi/MNDA/package/mnda")
`````
MNDA will also install TensorFlow and Keras for R, which need to be activated by installation of Miniconda. For this, according to the installation guidline of Keras ([here](https://cran.r-project.org/web/packages/keras/vignettes/index.html)):
`````{R}
library(keras)
install_keras()
`````
This is only required once for the installation.

## 2.2 Apply on simulated networks
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
We then feed ```graph_data``` to the MNDA pipleline specialized for a two-layer multiplex network (condition *"a"*), which is composed of two commands:
`````{R}
embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep = 50)
mnda_output = mnda_node_detection_2layer(embeddingSpaceList)
print(mnda_output$high_var_nodes_index)
`````
the ```mnda_embedding_2layer()``` function represents all the nodes in a common embedding space (step 1); and the ```mnda_node_detection_2layer()``` duncrion calculates the node-pair distances and assines a p-value to each node-pair (step 2 and 3). This process is repeated ```train.rep``` times to improve the robustness.

## 2.3. Usage Example 1: drug response  
*MNDA pipeline for condition "a"*

In this example, we construct gene coexpression netwerks (GCNs) for drug responders and non-responders. In this example, we load gene expression profile of cell lines of lung cancer as ```X``` and a binary vector of their response to the Tamoxifen drug as ```y```.
`````{R}
data = readRDS("Data/GCN2Layer_data_lung_tamoxifen_2000genes.rds")
X = data[[1]]
y = data[[2]]
`````
Next, we construc adjacency matrices of GCN for each condition,
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

## 2.4. Usage Example 2: application on individual specific networks


## Outreach
Any suggestions, collaboration and bug reports are welcome and appreciated.
Contact me via

Email: yousefi.bme@gmail.com,


[LinedIn](https://www.linkedin.com/in/behnam-yousefi-bme) or
[Twitter](https://twitter.com/behnam_bme)


