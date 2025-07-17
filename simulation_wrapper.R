library(msigdbr)
library(tidyverse)
library(FABLE)

get_geneSetMatrix = function(p, q){
  # msig_data = msigdbr(species = "Homo sapiens", category = "H")
  msig_data = msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'GO:MF')
  # select pathways
  gene_set = table(factor(msig_data$gene_symbol, levels = unique(msig_data$gene_symbol)),
                   msig_data$gs_name)
  gene_set = (gene_set > 0) * 1
  gene_set = gene_set[,sample(c(1:q), q, replace = F)]
  gene_set = gene_set[rowSums(gene_set) > 0,]
  gene_set = gene_set[sample(c(1:p), p, replace = F),]
  gene_set = gene_set[, colSums(gene_set) > 0]
  
  # N.B. it does not ensure q is the q selected at the beginning but q_selected<=q
  # I don't think is a problem, just put q bigger as input
  
  return(gene_set)
}


syntheticData = function(n, p, k, q, sigma_sq_0, sd_gamma, sd_psi, mseed, C_mispeficied=F){
  
  set.seed(mseed)
  C = get_geneSetMatrix(p, q)
  #C = unname(C)
  q = ncol(C)
  
  if(C_mispeficied==T){
    # set to zero the 20% of ones randomly selected
    i1 = which(C==1)
    idx0 = sample(x = i1, size = 0.2*length(i1), replace=F)
    C[idx0] = 0
    if(sum(colSums(C)==0)>0){
      C = C[, -which(colSums(C)==0)]
    }
    if(sum(rowSums(C)==0)>0){
      C = C[-which(rowSums(C)==0),]
    }
  }
  
  # compute P_C = C(C'C)^{-1}C' and Q_C = I_p - P_C
  tCCt = crossprod(C)
  s_cross_C = svd(tCCt)
  V_C = s_cross_C$u
  D_C = diag(as.vector(sqrt(s_cross_C$d)))
  D_C_inv = diag(as.vector(1/sqrt(s_cross_C$d)))
  U_C = C %*% V_C %*% D_C_inv
  P_C = tcrossprod(U_C)
  Q_C = diag(1, p, p) - P_C
  
  # factors
  M_0 = matrix(rnorm(n*k, 0, 1), ncol=k)
  # loadings
  Gamma_0 = matrix(rnorm(q*k, 0, sd_gamma), ncol=k)
  C_Gamma_0 = C %*% Gamma_0
  Psi_0 = Q_C %*% matrix(rnorm(p*k, 0, 0.2), ncol=k)
  Lambda_0 = C_Gamma_0 + Psi_0
  Lambda0_outer = tcrossprod(Lambda_0)
  
  # data
  Y = M_0 %*% t(Lambda_0) + sqrt(sigma_sq_0) * matrix(rnorm(n*p), nrow=n) 
  colnames(Y) = rownames(C)
  
  return(list("Y"=Y, "Lambda0_outer" = Lambda0_outer, "C"=C))
  
}


#----------------------------------# FABLE #----------------------------------#
# slight modification of code from https://github.com/shounakch/FABLE/tree/main
# to implement the FABLE methodology from
# Shounak Chattopadhyay, Anru R. Zhang, and David B. Dunson, (2024) 
# "Blessing of dimension in Bayesian inference on covariance matrices"
# arXiv preprint arXiv:2404.03805


PseudoPosteriorMean_2 <- function(Y,
                                  gamma0 = 1,
                                  delta0sq = 1,
                                  maxProp = 0.5) {
  
  tFABLEPostMean1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  
  kEst = RankEstimator(Y, 
                       U_Y,
                       V_Y,
                       svalsY,
                       kMax)
  
  FABLEHypPars = FABLEHyperParameters(Y,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst,
                                      gamma0,
                                      delta0sq)
  
  Part1 = FABLEHypPars$G
  Part2 = as.numeric(FABLEHypPars$gammaDeltasq / (FABLEHypPars$gamman - 2))
  CovEst = Part1 + diag(Part2)
  
  tFABLEPostMean2 = proc.time()
  tPostMean = (tFABLEPostMean2 - tFABLEPostMean1)[3]
  
  OutputList = list("FABLEPostMean" = CovEst,
                    "Lambda_outer" = Part1,
                    "FABLEHyperParameters" = FABLEHypPars,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "runTime" = tPostMean)
  
  return(OutputList)
  
}


PseudoPosteriorSampler_2 <- function(fit,
                                     Y,
                                     gamma0 = 1,
                                     delta0sq = 1,
                                     maxProp = 0.5,
                                     MC = 1000) {
  
  tFABLESample1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  #kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  kEst = fit$estRank
  
  FABLEHypPars = FABLEHyperParameters(Y,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst,
                                      gamma0,
                                      delta0sq)
  
  CovCorrectMatrix = cov_correct_matrix(FABLEHypPars$SigmaSqEstimate, 
                                        FABLEHypPars$G)
  
  varInflation = (sum(CovCorrectMatrix) / (p*(p+1)/2))^2
  
  FABLESamples = FABLESampler(Y, 
                              gamma0, 
                              delta0sq, 
                              MC,
                              U_Y,
                              V_Y,
                              svalsY,
                              kEst,
                              FABLEHypPars$tauSqEstimate,
                              FABLEHypPars$gammaDeltasq,
                              FABLEHypPars$G0,
                              varInflation)
  
  tFABLESample2 = proc.time()
  tSample = (tFABLESample2 - tFABLESample1)[3]
  
  OutputList = list("CCFABLESamples" = FABLESamples,
                    "FABLEHyperParameters" = FABLEHypPars,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "varInflation" = varInflation,
                    "runTime" = tSample)
  
  return(OutputList)
  
}