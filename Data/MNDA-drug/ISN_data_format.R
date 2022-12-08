### Convert ISNs to the standard MNDA file format ###
# This code is data-specific and can be modified for the specific use
# The current script is written for the project: MNDA-drug

rm(list = ls())
library(data.table)
library(stringr)

setwd("~/Desktop/R_Root/MNDA")

### 1) Convert ISNs to the standard MNDA file format ###
data = data.frame(fread(file="Data/MNDA-drug/CD_TNF_w14_clinical/ISN_Resulting_net_notNULL_NONODE_MAGMACONF.txt",sep = " " ))
Nodelist = t(sapply(data[,1], function(x) strsplit(x,"_")[[1]]))
IndvGraphWeights = data[,-1]
IndvGraphWeights = abs(IndvGraphWeights)
IndvGraphWeights = (IndvGraphWeights - min(IndvGraphWeights)) / (max(IndvGraphWeights) - min(IndvGraphWeights))

outcome = read.delim("Data/MNDA-drug/CD_TNF_w14_clinical/w14_n188_metadata.txt")

rownames(outcome) = str_replace(outcome$FC.nummer, "-", ".")
indv = intersect(rownames(outcome), colnames(IndvGraphWeights))

outcome = outcome[indv, ]
IndvGraphWeights = IndvGraphWeights[,indv]

outcome$response = outcome$Clinical_outcome_combined       # for CD_TNF_w14_clinical
# outcome$response = outcome$Endoscopic_outcome_combined     # for CD_TNF_w24_endoscopic

data_new = cbind(Nodelist, IndvGraphWeights)
outcome_new = outcome$response

# saveRDS(data_new, file = "Data/MNDA-drug/CD_TNF_w14_ISN.rds")
# saveRDS(outcome_new, file = "Data/MNDA-drug/CD_TNF_w14_outcome.rds")

## 2) Convert ISNs to several global networks based on the a given variable ###
data = data.frame(fread(file="Data/MNDA-drug/CD_TNF_w14_clinical/ISN_Resulting_net_notNULL_NONODE_MAGMACONF.txt",sep = " " ))
Nodelist = t(sapply(data[,1], function(x) strsplit(x,"_")[[1]]))
IndvGraphWeights = data[,-1]

outcome = read.delim("Data/MNDA-drug/CD_TNF_w14_clinical/w14_n188_metadata.txt")

rownames(outcome) = str_replace(outcome$FC.nummer, "-", ".")
indv = intersect(rownames(outcome), colnames(IndvGraphWeights))

outcome = outcome[indv, ]
IndvGraphWeights = IndvGraphWeights[,indv]

outcome$response = outcome$Clinical_outcome_combined       # for CD_TNF_w14_clinical
# outcome$response = outcome$Endoscopic_outcome_combined     # for CD_TNF_w24_endoscopic

global_net_responder = apply(IndvGraphWeights[outcome$response == 1], 1, mean)
global_net_nonresponder = apply(IndvGraphWeights[outcome$response == 0], 1, mean)

global_net_weights = cbind(global_net_responder, global_net_nonresponder)
global_net_weights = abs(global_net_weights)
global_net_weights = (global_net_weights - min(global_net_weights)) / (max(global_net_weights) - min(global_net_weights))


data_new = cbind(Nodelist, global_net_weights)

saveRDS(data_new, file = "Data/MNDA-drug/CD_TNF_w14_Global.rds")

## Node labels
Tab = read.delim("Data/MNDA-drug/CD_TNF_w14_clinical/Sequences_nodes_complete.txt", sep = " ")
rownames(Tab) = Tab$Node_list

Names1 = Tab[,"Family"]
Names2 = Tab[,"Genus"]
Names1[is.na(Names1)] = "unclass. family"
Names2[is.na(Names2)] = "genus"
Names = paste(Names1, Names2)

Tab$Labels = Names
Tab = Tab[names(V(graph)),]

high_var_nodes = c(130, 124, 131, 136, 104, 115, 135,  99, 134, 112, 40)



