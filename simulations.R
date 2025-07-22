library(PLIER)

rm(list=ls())

source('helper_functions.R')
source("simulation_wrapper.R")
source("FACTOR_ANALYSIS/FACTOR_CODE_update.R")

# set parameters
Nsim = 25
param = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 0.5,
             sd_gamma = 2, sd_psi = 0.5)
# param = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 0.5,
#              sd_gamma = 1, sd_psi = 0.3)


set.seed(463)
seeds_g = sample.int(9000, Nsim)

kfitBASIL = err_normBASIL = timeBASIL = rep(0, Nsim)
# kfitFABLE = err_normFABLE  = timeFABLE = rep(0, Nsim)
kfitROTATE = err_normROTATE  = timeROTATE = rep(0, Nsim)
kfitPLIER = err_normPLIER  = timePLIER = rep(0, Nsim)
data_all = vector("list", Nsim)

for(s in 1:Nsim){
  #simulate data
  datas = syntheticData(n = param$n, p = param$p, k = param$k, q = param$q, 
                        sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                        sd_psi = param$sd_psi, mseed = seeds_g[s], C_mispeficied = F)
  data_all[[s]] = datas
  Ys = datas$Y
  Cs = datas$C
  Lambda0_outer = datas$Lambda0_outer
  
  # compute BASIL
  ptmB = proc.time()
  est_kBASIL = estimate_latent_dimension(Ys, k_max = 50)
  fitBASIL = compute_point_estimates(Ys, Cs, k = est_kBASIL$k_hat)
  etmB = proc.time() - ptmB
  timeBASIL[s] = etmB[1] + etmB[2]
  Lambda_BASIL = fitBASIL$Lambda_C + fitBASIL$Lambda_N
  err_normBASIL[s] = norm(tcrossprod(Lambda_BASIL) - Lambda0_outer, type='F')/
    norm(Lambda0_outer, type='F')
  kfitBASIL[s] = est_kBASIL$k_hat

  # compute rotate

  K = 30
  startB = matrix(rnorm(param$p*K),param$p,K)
  start = list(B = startB, sigma = rep(1,K), theta = rep(0.5,K))
  ptmR = proc.time()
  fitROTATE = FACTOR_ROTATE(Y = Ys, lambda0 = 5, lambda1= 0.001, start = start, 
                            K = K, epsilon = 0.05, alpha = 1/param$p, PX = TRUE,
                            approximate = TRUE, stop = 100, varimax = TRUE, plot = FALSE)
  etmR = proc.time() - ptmR
  timeROTATE[s] = etmR[1] + etmR[2]
  Lambda_ROTATE = fitROTATE$B
  err_normROTATE[s] = norm(tcrossprod(Lambda_ROTATE) - Lambda0_outer, type='F')/
    norm(Lambda0_outer, type='F')
  
  # # compute FABLE
  # ptmF = proc.time()
  # fitFABLE = PseudoPosteriorMean_2(Ys)
  # etmF = proc.time() - ptmF
  # timeFABLE[s] = etmF[1] + etmF[2]
  # kfitFABLE[s] = fitFABLE$estRank
  # err_normFABLE[s] = norm(fitFABLE$Lambda_outer - Lambda0_outer, type='F')/
  #   norm(Lambda0_outer, type='F')
  # 
  # # compute PLIER
  # ptmP = proc.time()
  # fitPLIER = PLIER(t(Ys), Cs, scale = F, minGenes = 1, doCrossval = T)
  # etmP = proc.time() - ptmP
  # timePLIER[s] = etmP[1] + etmP[2]
  # kfitPLIER[s] = nrow(fitPLIER$B)
  # covPLIER = (fitPLIER$Z %*% (fitPLIER$B %*% t(fitPLIER$B)) %*% t(fitPLIER$Z))/param$n
  # err_normPLIER[s] = norm(covPLIER - Lambda0_outer, type='F')/
  #   norm(Lambda0_outer, type='F')
  
}

summary(err_normROTATE)
summary(err_normBASIL)

# c(median(kfitBASIL), IQR(kfitBASIL))
# c(median(kfitFABLE), IQR(kfitFABLE))
# c(median(err_normBASIL), IQR(err_normBASIL))
# c(median(err_normROTATE), IQR(err_normBASIL))
# c(median(err_normFABLE), IQR(err_normFABLE))
# c(median(timeBASIL), IQR(timeBASIL))
# c(median(timeFABLE), IQR(timeFABLE))
# 
# c(median(kfitPLIER), IQR(kfitPLIER))
# c(median(err_normPLIER), IQR(err_normPLIER))
# c(median(timePLIER), IQR(timePLIER))
# 
# SimRes1 = cbind.data.frame(c(err_normBASIL, err_normFABLE, err_normPLIER), 
#                           c(timeBASIL, timeFABLE, timePLIER),
#                           c(kfitBASIL, kfitFABLE, kfitPLIER),
#                           c(rep("BASIL",Nsim), rep("FABLE",Nsim),  rep("PLIER",Nsim)),
#                           rep("3000", Nsim*3))
# colnames(SimRes1) = c("err_norm", "time", "k_est", "model", "p")
# 
# # save(SimRes, file="simResults.Rdata")

