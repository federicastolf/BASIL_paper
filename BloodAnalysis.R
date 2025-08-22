library(PLIER)
library(tidyverse)
library(mvtnorm)
rm(list=ls())

#---------------------------# human whole blood data #------------------------#

data("dataWholeBlood") 
data("bloodCellMarkersIRISDMAP")
data("canonicalPathways")
allPaths = combinePaths(bloodCellMarkersIRISDMAP, canonicalPathways)
cm=intersect(rownames(dataWholeBlood), rownames(allPaths))
allPaths=allPaths[cm,]
dataWholeBlood=dataWholeBlood[cm,]
Y = t(dataWholeBlood)
Y = scale(Y)
n = nrow(Y)

#----# fit BASIL #-----#
source('helper_functions.R')
est_k = estimate_latent_dimension(Y, k_max=30)
fitBL = compute_point_estimates(Y,  allPaths, k=est_k$k_hat)
Lambda_hat = fitBL$Lambda_C + fitBL$Lambda_N
CovBasil = Lambda_hat %*% t(Lambda_hat) + fitBL$sigma_sq*diag(ncol(Y))
corrBasil = cov2cor(CovBasil)


#----------# out of sample prediction #---------#

set.seed(423)
train_set = sample(1:n, as.integer(0.8*n))
Y_train = Y[train_set,]
Y_test = Y[-train_set,]

#-----# fit BASIL
est_k = estimate_latent_dimension(Y_train, k_max=25)
train_fit = compute_point_estimates(Y_train, allPaths, k=est_k$k_hat)
Lambda_hat = train_fit$Lambda_C + train_fit$Lambda_N
CovBasil = Lambda_hat %*% t(Lambda_hat) + train_fit$sigma_sq*diag(ncol(Y))

#----# fit PLIER 
library(matrixcalc)
PLresult = PLIER(t(Y_train), allPaths)
covPLIER = PLresult$Z %*% (cov(t(PLresult$B))) %*% t(PLresult$Z)
#covPLIER = (PLresult$Z %*% (PLresult$B %*% t(PLresult$B)) %*% t(PLresult$Z))/nrow(Y_train)
is.positive.definite(covPLIER)
isSymmetric.matrix(covPLIER)
c = 0.5*(covPLIER+t(covPLIER)) + diag(rep(1e-7, ncol(covPLIER)))
# isSymmetric.matrix(c)

basil_ll_00s = dmvnorm(Y_test, sigma = CovBasil, log=T)
plier_ll_00s = dmvnorm(Y_test, sigma = c, log=T)
#plier_ll_00s = dmvnorm(Y_test, sigma = covPLIER, log=T)
t.test(basil_ll_00s, plier_ll_00s , alternative='greater', paired=TRUE)

