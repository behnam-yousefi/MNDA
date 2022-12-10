# Multiplex Network Differential Analysis (MNDA) <img width="100" alt="image" src="https://github.com/behnam-yousefi/MNDA/blob/master/Figures/MNDA_logo.png?raw=true">

Interactions between different biological elements are crucial for the function of biological systems and are usually represented as networks. The dysregulation of these networks can be associated with different clinical conditions such as diseases and response to treatments. In this work,  we propose an R package, called multiplex network differential analysis (MNDA) to quantify and test the variations in the local neighborhood of nodes between the two given conditions (e.g. case-control). We further show examples of finding important subnetworks in gene co-expression networks for response to treatment as well as an use case for individual specific networks (ISNs).

The core of the MNDA tool consists of three steps (Figure 1):
1. representing the nodes of all networks layers into a common embedding space (using EDNN);
2. calculate the distance between the nodes corresponding to the same element (e.g. gene);
3. detect the nodes whose neighborhood varies significantly based on statistical testing (using permuted graphs).

<img width="800" alt="image" src="https://github.com/behnam-yousefi/MNDA/blob/master/Figures/Figure_1.png?raw=true">

**Figure 2.** The schematic representation of the MNDA workflow. All the nodes of all the layers along with the permuted networks are represented into a common embedding space. The distances between all the pairs  of the permitted network are used to construct a null probability distribution fiction (PDF), based on which statistical testing is performed to detect the nodes whose neighborhood significantly changes.


The current MNDA pipeline is desined for two conditions:
a. two-layer network case corresponding to two (paired/unpaired) conditions (e.g. healthy-disease);
b. multi-layer network case (e.g. ISNs) with two matched groups (e.g. before treatment-after treatment).

The future features to be considered:
c. multi-layer network case (e.g. ISNs) with two unmatched groups (e.g. healthy-disease);
d. multi-layer network case each layer with different condition (e.g. temporal networks).


**Note:** the first step (EDNN) is the same for all the cases but they are different in step 2 and 3.

**Note:** the network structures for different layers shoul be the same but different in weights.

Reference: [to be added]

## Installation
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

## Simulated Network
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
embeddingSpaceList = mnda_embedding_2layer(graph_data, train.rep = 50, walk.rep = 100, random.walk = TRUE, null.perm = TRUE)

Results = mnda_node_detection_2layer(embeddingSpaceList, p.adjust.method = "none")
print(Results$high_var_nodes_index)
`````
the ```mnda_embedding_2layer()``` function represents all the nodes in a common embedding space (step 1).
the ```mnda_node_detection_2layer()``` duncrion calculates the node-pair distances and assines a p-value to each node-pair (step 2 and 3).

## Usage Example 1: drug response  

## Usage Example 2: application on individual specific networks
