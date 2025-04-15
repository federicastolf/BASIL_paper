library(PLIER)
library(msigdbr)

rm(list=ls())

# Load Data
# data2 = read.csv("SICK-AIM3/Data/logCPM.csv")
dataSICK = data2[,-1]
rownames(dataSICK) <- data2[,1]

# create prior knowledge matrix
GOSet = msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'GO:MF')
geneCats = unique(GOSet$gs_name)
C = matrix(nrow = nrow(dataSICK), ncol = length(geneCats))
colnames(C) = geneCats
row.names(C) = row.names(dataSICK)

for(j in 1:ncol(C)){
  C[,j] = as.numeric(row.names(C) %in% GOSet$ensembl_gene[GOSet$gs_name==colnames(C)[j]])
}

cc = colSums(C)
idc0 = which(cc==0)
C = C[,-idc0]
rc = rowSums(C)
idr0 = which(rc==0)
C = C[-idr0,]
dataSICK = dataSICK[-idr0 ,]
dataSICK= t(dataSICK)
dataSICK = scale(dataSICK)

save(C, dataSICK, file="DataSICK.Rdata")


