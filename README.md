# Multiplex Network Differential Analysis (MNDA)

logo

Interactions between different biological elements are crucial for the function of biological systems and are usually represented as networks. The dysregulation of these networks can be associated with different clinical conditions such as diseases and response to treatments. In this work,  we propose an R package, called multiplex network differential analysis (MNDA) to quantify and test the variations in the local neighborhood of nodes between the two given conditions (e.g. case-control). We further show examples of finding important subnetworks in gene co-expression networks for response to treatment as well as an use case for individual specific networks (ISNs).

The core of the MNDA tool consists of three steps:
1. representing the nodes of all networks layers into a common embedding space (using EDNN);
2. calculate the distance between the nodes corresponding to the same element (e.g. gene);
3. detect the nodes whose neighborhood varies significantly based on statistical testing (using permuted graphs).

The current MNDA pipeline is desined for:
* two-layer network case corresponding to two (paired/unpaired) conditions (e.g. healthy-disease);
* multi-layer network case (e.g. ISNs) with two matched groups (e.g. before treatment-after treatment).

The future features to be considered:
* multi-layer network case (e.g. ISNs) with two unmatched groups (e.g. healthy-disease);
* multi-layer network case each layer with different condition (e.g. temporal networks).


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

## Usage Example 1: drug response  

## Usage Example 2: application on individual specific networks
