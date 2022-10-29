# Convert ISNs to the standard format

data = data.frame(fread(file="Data/Resulting_net_notNULL_MAGMACONF6M.txt",sep = " " ))
Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
MianGraph1 = graph(Nodelist, directed = FALSE)
IndvGraphWeights1 = data[,-1]
IndvGraphWeights1 = abs(IndvGraphWeights1)
IndvGraphWeights1 = (IndvGraphWeights1 - min(IndvGraphWeights1)) / (max(IndvGraphWeights1) - min(IndvGraphWeights1))

data = data.frame(fread(file="Data/Resulting_net_notNULL_MAGMACONF9M.txt",sep = " " ))
Nodelist = sapply(data$X, function(x) strsplit(x,"_")[[1]])
MianGraph2 = graph(Nodelist, directed = FALSE)
IndvGraphWeights2 = data[,-1]
IndvGraphWeights2 = abs(IndvGraphWeights2)
IndvGraphWeights2 = (IndvGraphWeights2 - min(IndvGraphWeights2)) / (max(IndvGraphWeights2) - min(IndvGraphWeights2))

LuckiMap1 = read.delim("Data/LuckiMap_6M.txt")
LuckiMap2 = read.delim("Data/LuckiMap_9M.txt")
Children =  merge(LuckiMap1, LuckiMap2, by = "Child")

IndvGraphWeights1 = IndvGraphWeights1[,Children$Sample.name.x]
IndvGraphWeights2 = IndvGraphWeights2[,Children$Sample.name.y]