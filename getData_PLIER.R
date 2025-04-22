
#install_github("wgmao/PLIER")
library(PLIER)

rm(list=ls())


#---------------------------# human whole blood data #------------------------#

data("dataWholeBlood") # omics data
data("bloodCellMarkersIRISDMAP")
data("canonicalPathways")

allPaths = combinePaths(bloodCellMarkersIRISDMAP, canonicalPathways)
# then they do the intersection between the two dataset for finding common genes
# (They also consider a minimum number of genes a pathway must have to be considered)

cm=intersect(rownames(dataWholeBlood), rownames(allPaths))
allPaths=allPaths[cm,]
dataWholeBlood=dataWholeBlood[cm,]

dataWholeBlood = t(dataWholeBlood)
summary(c(dataWholeBlood))
# not standardize (mean 0 and sd1) -> maybe we should do it?
# in PLIER is the default, so seems like standardize make sense
dataWholeBloods = scale(dataWholeBlood)
# eventually we obtain C 5892x606 and Y 36x5892
# save(dataWholeBloods, allPaths, file="dataBlood.Rdata")

#-------------------# vaccination dataset #----------------#
rm(list = ls())

data("bloodCellMarkersIRISDMAP")
data("svmMarkers")
data("canonicalPathways")
data("vacData")

allPaths = combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, canonicalPathways)

cm=intersect(rownames(vacData), rownames(allPaths))
allPaths=allPaths[cm,]
vacData=vacData[cm,]
# we obtain C 5359x628 and Y 163x5359

vacData = t(vacData)
vacData = scale(vacData)

# save(vacData, allPaths, file="dataVac.Rdata")

#-------------------# single-cell RNA seq from Usokin et al. #----------------#
rm(list = ls())

data("chemgenPathways")
data("canonicalPathways")
data("dataUsoskin")
data("human2Mouse")
data("cellsUsoskin")

# Create the combined pathway matrix
usoskinPath = combinePaths(canonicalPathways, chemgenPathways)
# Map the pathway matrix to mouse names
usoskinPath = mapPathway(usoskinPath, human2Mouse)

cm = commonRows(dataUsoskin, usoskinPath)

dataUsoskin = dataUsoskin[cm,]
usoskinPath = usoskinPath[cm,] # C: 8710 x 3940
dataUsoskin = rowNorm(dataUsoskin)
dataUsoskin = t(dataUsoskin) # Y:622 x 8710

# save(dataUsoskin, usoskinPath, file="Usoskindata.Rdata")
