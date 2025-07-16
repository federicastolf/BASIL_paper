library(FABLE)
library(PLIER)

rm(list=ls())

source('helper_functions.R')
source("simulation_wrapper.R")

# set parameters
Nsim = 25
param = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 0.5,
             sd_gamma = 1, sd_psi = 0.3)
# param = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 0.5, 
#              sd_gamma = 1, sd_psi = 0.3)
  
  
set.seed(463)
seeds_g = sample.int(9000, Nsim)

kfitBASIL = err_normBASIL = timeBASIL = rep(0, Nsim)
kfitFABLE = err_normFABLE  = timeFABLE = rep(0, Nsim)
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

  # compute FABLE
  ptmF = proc.time()
  fitFABLE = FABLE::PseudoPosteriorMean(Ys)
  etmF = proc.time() - ptmF
  timeFABLE[s] = etmF[1] + etmF[2]
  kfitFABLE[s] = fitFABLE$estRank
  err_normFABLE[s] = norm(fitFABLE$FABLEPostMean - Lambda0_outer, type='F')/
    norm(Lambda0_outer, type='F')
  
  # compute PLIER
  ptmP = proc.time()
  fitPLIER = PLIER(t(Ys), Cs, scale = F, minGenes = 1, doCrossval = T)
  etmP = proc.time() - ptmP
  timePLIER[s] = etmP[1] + etmP[2]
  kfitPLIER[s] = nrow(fitPLIER$B)
  covPLIER = fitPLIER$Z %*% (fitPLIER$B %*% t(fitPLIER$B)) %*% t(fitPLIER$Z)
  err_normPLIER[s] = norm(covPLIER - Lambda0_outer, type='F')/
    norm(Lambda0_outer, type='F')
  
}

c(median(kfitBASIL), IQR(kfitBASIL))
c(median(kfitFABLE), IQR(kfitFABLE))
c(median(err_normBASIL), IQR(err_normBASIL))
c(median(err_normFABLE), IQR(err_normFABLE))
c(median(timeBASIL), IQR(timeBASIL))
c(median(timeFABLE), IQR(timeFABLE))

c(median(kfitPLIER), IQR(kfitPLIER))
c(median(err_normPLIER), IQR(err_normPLIER))
c(median(timePLIER), IQR(timePLIER))


