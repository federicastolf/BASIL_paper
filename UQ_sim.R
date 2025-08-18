rm(list=ls())

source('helper_functions.R')
source("simulation_wrapper.R")

# set parameters
Nsim = 10
# high biological signal
# param = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7,
#              sd_psi = 0.1)

#low biological signal
param = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.4,
             sd_psi = 0.7)

set.seed(463)
seeds_g = sample.int(9000, Nsim)
ccCoverage = rep(0,Nsim)

for(s in 1:Nsim){
  #simulate data
  datas = syntheticData(n = param$n, p = param$p, k = param$k, q = param$q, 
                        sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                        sd_psi = param$sd_psi, mseed = seeds_g[s])
  Ys = datas$Y
  Cs = datas$C
  Lambda0_outer = datas$Lambda0_outer
  
  # compute BASIL
  fitBASIL = compute_point_estimates(Ys, Cs, k = param$k)
  subsample_index = 1:200
  alpha = 0.05
  
  ## sample loadings given M = M_hat
  params_posterior_samples = compute_posterior_samples_cc(
    Ys, fitBASIL$Lambda_C, fitBASIL$Lambda_N, fitBASIL$tau_C, fitBASIL$tau_N, 
    fitBASIL$sigma_sq, fitBASIL$P_C)
  
  Lambda_outer_posterior_samples = sample_Lambda_outer(
    params_posterior_samples$Lambda_samples[subsample_index,,])
  
  Lambda_outer_qs = apply(Lambda_outer_posterior_samples, c(1,2), 
                          function(x)(quantile(x, probs=c(alpha/2, 1-alpha/2))))
  cov = mean((Lambda_outer_qs[1,,]<Lambda0_outer[subsample_index, subsample_index]) &
         (Lambda_outer_qs[2,,]>Lambda0_outer[subsample_index, subsample_index]))
  ccCoverage[s] = cov
}

summary(ccCoverage)
