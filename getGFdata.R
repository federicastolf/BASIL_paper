
library(tidyverse)
library(limma)
# library(edgeR)
library(openxlsx)
library(readxl)
library(data.table)
library(org.Hs.eg.db)
library(msigdbr)
library(AnnotationDbi)


#--------------------------------# Load Data #---------------------------------#
data <- read.table(file = 'GSE211567_normData_discovery_2021MAR24.txt', 
                   sep = '\t', row.names = 1, header = TRUE)
data[1:5, 1:5]

#----------------------# Construct Prior Knowledge Matrix #--------------------#

# map between Ensembl IDs and RefSeq IDs
map <- AnnotationDbi::select(org.Hs.eg.db, key = row.names(data), columns = c("ENSEMBL"),
    keytype = "REFSEQ")
map <- map[is.na(map$ENSEMBL) == FALSE, ]
length(unique(map$ENSEMBL))

# get gene set information
GOSet <- msigdbr(species = 'Homo sapiens', collection = 'C5', subcollection = 'GO:MF')
head(GOSet)
GOSet <- merge(GOSet, map, by.x = 'ensembl_gene', by.y = 'ENSEMBL', sort = FALSE)

# construct prior knowledge matrix
lcpm <- data
geneCats <- unique(GOSet$gs_name)
geneSetMat <- matrix(nrow = nrow(lcpm), ncol = length(geneCats))
colnames(geneSetMat) <- geneCats
row.names(geneSetMat) <- row.names(lcpm)

for(j in 1:ncol(geneSetMat)){
  geneSetMat[ ,j] <- as.numeric(row.names(geneSetMat) %in% 
                                  GOSet$REFSEQ[GOSet$gs_name == colnames(geneSetMat)[j]])
  if((j %% 100)==0){print(paste0('Iteration ', j, ' complete'))}
}

# filtered out pathways with fewer than 10 genes and genes that were not 
# annotated to any gene set
sc = colSums(geneSetMat)
pr1 = which(sc<11)
geneSetMat = geneSetMat[,-pr1]
sr = rowSums(geneSetMat)
gr1 = which(sr==0)
geneSetMat = geneSetMat[-gr1,]
data = data[-gr1,]

# save(data, geneSetMat, file="data_Gfever.Rdata")
