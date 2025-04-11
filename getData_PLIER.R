
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

