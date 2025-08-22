library(PLIER)
rm(list = ls())

# load data and functions
data("bloodCellMarkersIRISDMAP")
data("svmMarkers")
data("canonicalPathways")
data("vacData")
source('helper_functions.R')

#----# obtain data 
allPaths = combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, canonicalPathways)
cm=intersect(rownames(vacData), rownames(allPaths))
allPaths=allPaths[cm,]
vacData=vacData[cm,]
# we obtain C 5359x628 and Y 163x5359

vacData = t(vacData)
vacData = scale(vacData)
n = nrow(vacData)

#---------# fit BASIL #------------#
est_k = estimate_latent_dimension(vacData , k_max = 50)
vacData_fit = compute_point_estimates(vacData, allPaths, k = est_k$k_hat)
vacData_fit$tau_C/vacData_fit$tau_N

Lambda_BASIL = vacData_fit$Lambda_C + vacData_fit$Lambda_N
cov_vacData = tcrossprod(Lambda_BASIL) +vacData_fit$sigma_sq*diag(ncol(vacData))
corrBasil = cov2cor(cov_vacData)
summary(c(corrBasil))


#----------# out of sample prediction #---------#

set.seed(423)
train_set = sample(1:n, as.integer(0.8*n))
Y_train = vacData[train_set,]
Y_test = vacData[-train_set,]

#-----# fit BASIL
est_k = estimate_latent_dimension(Y_train, k_max=25)
train_fit = compute_point_estimates(Y_train, allPaths, k=est_k$k_hat)
Lambda_hat = train_fit$Lambda_C + train_fit$Lambda_N
CovBasil = Lambda_hat %*% t(Lambda_hat) + train_fit$sigma_sq*diag(ncol(vacData))

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

