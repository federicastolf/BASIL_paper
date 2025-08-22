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
est_kWB = estimate_latent_dimension(Y, k_max=30)
fitWB = compute_point_estimates(Y,  allPaths, k=est_kWB$k_hat)
fitWB$tau_C/fitWB$tau_N

Lambda_hatWB = fitWB$Lambda_C + fitWB$Lambda_N
CovBasilWB = Lambda_hatWB %*% t(Lambda_hatWB) + fitWB$sigma_sq*diag(ncol(Y))
corrBasilWb = cov2cor(CovBasilWB)
summary(c(corrBasilWb))

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
PLresult = PLIER(t(Y_train), allPaths, scale = F)
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
 
