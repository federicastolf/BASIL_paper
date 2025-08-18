library(PLIER)

rm(list=ls())

source('helper_functions.R')
source("simulation_wrapper.R")

# set parameters
Nsim = 25
# high biological signal
# param = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7,
#              sd_psi = 0.1)

#low biological signal
param = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.4,
             sd_psi = 0.7)


set.seed(463)
seeds_g = sample.int(9000, Nsim)

kfitBASIL = err_normBASIL = timeBASIL = rep(0, Nsim)
kfitROTATE = err_normROTATE  = timeROTATE = rep(0, Nsim)
kfitPLIER = err_normPLIER  = timePLIER = rep(0, Nsim)
data_all = vector("list", Nsim)

for(s in 1:Nsim){
  #simulate data
  datas = syntheticData(n = param$n, p = param$p, k = param$k, q = param$q, 
                        sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                        sd_psi = param$sd_psi, mseed = seeds_g[s])
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

  # compute ROTATE
  K = param$k
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
  
  # compute PLIER
  ptmP = proc.time()
  fitPLIER = PLIER(t(Ys), Cs, scale = F, minGenes = 1, doCrossval = T)
  etmP = proc.time() - ptmP
  timePLIER[s] = etmP[1] + etmP[2]
  kfitPLIER[s] = nrow(fitPLIER$B)
  covPLIER = (fitPLIER$Z %*% (fitPLIER$B %*% t(fitPLIER$B)) %*% t(fitPLIER$Z))/param$n
  err_normPLIER[s] = norm(covPLIER - Lambda0_outer, type='F')/
    norm(Lambda0_outer, type='F')

}

summary(err_normROTATE)
summary(err_normBASIL)
summary(err_normPLIER)


# SimResLow1 = cbind.data.frame(c(err_normBASIL, err_normROTATE, err_normPLIER),
#                           c(timeBASIL, timeROTATE, timePLIER),
#                           c(kfitBASIL, kfitROTATE, kfitPLIER),
#                           c(rep("BASIL",Nsim), rep("ROTATE",Nsim),  rep("PLIER",Nsim)),
#                           rep("3000", Nsim*3))
# colnames(SimResLow1) = c("err_norm", "time", "k_est", "model", "p")
# # 
# 
# load("simResultsLow.Rdata")
# SimResLow = rbind.data.frame(SimResLow, SimResLow1)
# save(SimResLow, file="simResultsLow.Rdata")
# 

# B1 = SimResLow %>% filter(p=="3000") %>% filter(model=="PLIER") 
# c(median(B1$err_norm), IQR(B1$err_norm))
# c(median(B1$k_est), IQR(B1$k_est))
# c(median(B1$time), IQR(B1$time))

