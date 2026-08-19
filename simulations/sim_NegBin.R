library(PLIER)
library(BASIL)
library(parallel)

rm(list = ls())

source("functs/FACTOR_CODE_update.R")   
source("functs/helper.R")              

# eror helper
rel_err = function(A, B) norm(A - B, type = 'F') / norm(B, type = 'F')

cor_from_outer = function(Lam, sigma_sq) {
  d = diag(Lam) + sigma_sq
  d[d <= 0] = .Machine$double.eps
  R = Lam * tcrossprod(1/sqrt(d))
  diag(R) = 1
  R
}


#------------# accuracy simulation helper - BASIL, ROTATE, PLIER #-------------#

run_simulation_studyNB = function(param, scenario_name, Nsim, seed, ncores = 6) {
  
  set.seed(seed)
  seeds_g = sample.int(9000, Nsim)
  
  one_sim = function(s) {
    datas = syntheticDataNB(n = param$n, p = param$p, k = param$k, q = param$q,
                            sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                            sd_psi = param$sd_psi, mu_0 = param$mu_0,
                            size_nb = param$size_nb, mseed = seeds_g[s])
    zero_frac = mean(datas$Y_counts == 0)
    
    # log1p counts, column centred
    Ys = matrix(scale(datas$Y), nrow = param$n)
    colnames(Ys) = colnames(datas$Y)
    Cs = datas$C
    Lambda0_outer = datas$Lambda0_outers
    
    R0 = cor_from_outer(Lambda0_outer, param$sigma_sq_0)
    
    ## BASIL
    ptmB = proc.time()
    est_kBASIL = estimate_latent_dimension(Ys, k_max = 50)
    fitBASIL = BASIL_point_estimates(Ys, Cs, k = est_kBASIL$k_hat)
    etmB = proc.time() - ptmB
    timeBASIL = etmB[1] + etmB[2]
    Lambda_BASIL = tcrossprod(fitBASIL$Lambda_C + fitBASIL$Lambda_N)
    corr_BASIL = cor_from_outer(Lambda_BASIL, fitBASIL$sigma_sq)
    kfitBASIL    = est_kBASIL$k_hat
    err_corBASIL = rel_err(corr_BASIL, R0)
    
    ## ROTATE
    K = param$k
    startB = matrix(rnorm(param$p * K), param$p, K)
    start = list(B = startB, sigma = rep(1, K), theta = rep(0.5, K))
    ptmR = proc.time()
    fitROTATE = FACTOR_ROTATE(Y = Ys, lambda0 = 5, lambda1 = 0.001, start = start,
                              K = K, epsilon = 0.05, alpha = 1 / param$p, PX = TRUE,
                              approximate = TRUE, stop = 100, varimax = TRUE, plot = FALSE)
    etmR = proc.time() - ptmR
    timeROTATE = etmR[1] + etmR[2]
    Lambda_ROTATE = tcrossprod(fitROTATE$B)
    corr_ROTATE = cor_from_outer(Lambda_ROTATE, mean(fitROTATE$sigma^2))
    kfitROTATE    = ncol(fitROTATE$B)
    err_corROTATE = rel_err(corr_ROTATE, R0)
    
    ## PLIER
    ptmP = proc.time()
    fitPLIER = PLIER(t(Ys), Cs, scale = FALSE, minGenes = 1, doCrossval = TRUE)
    etmP = proc.time() - ptmP
    timePLIER = etmP[1] + etmP[2]
    kfitPLIER = nrow(fitPLIER$B)
    covPLIER = (fitPLIER$Z %*% tcrossprod(fitPLIER$B) %*% t(fitPLIER$Z)) / param$n
    sigma_sqPLIER = mean((Ys - t(fitPLIER$B) %*% t(fitPLIER$Z))^2)
    corr_PLIER = cor_from_outer(covPLIER, sigma_sqPLIER)
    err_corPLIER = rel_err(corr_PLIER, R0)
    
    data.frame(
      err_cor   = c(err_corBASIL, err_corROTATE, err_corPLIER),
      time      = c(timeBASIL,    timeROTATE,    timePLIER),
      k_est     = c(kfitBASIL,    kfitROTATE,    kfitPLIER),
      model     = c("BASIL", "ROTATE", "PLIER"),
      zero_frac = zero_frac,
      sim       = s,
      stringsAsFactors = FALSE)
  }
  
  res_list = mclapply(1:Nsim, one_sim, mc.cores = ncores, mc.preschedule = FALSE)
  
  df_combined = do.call(rbind, res_list[ok])
  df_combined$scenario = scenario_name
  df_combined$p = as.character(param$p)
  rownames(df_combined) = NULL
  
  cat(sprintf("%s: completed %d/%d simulations\n", scenario_name, sum(ok), Nsim))
  return(df_combined)
}


#------------------# coverage simulations BASIL helper #-----------------------#


run_coverage_simulationNB = function(param, scenario_name, subsample_index, alpha,
                                     Nsim, seed, ncores = 6) {
  
  set.seed(seed)
  seeds_g = sample.int(9000, Nsim)

  one_sim = function(s) {
    datas = syntheticDataNB(n = param$n, p = param$p, k = param$k, q = param$q,
                            sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                            sd_psi = param$sd_psi, mu_0 = param$mu_0,
                            size_nb = param$size_nb, mseed = seeds_g[s])
    zero_frac = mean(datas$Y_counts == 0)
    
    # log1p counts, column centred
    Ys = matrix(scale(datas$Y), nrow = param$n)
    colnames(Ys) = colnames(datas$Y)
    Cs = datas$C
    
    # target: true correlation restricted to the subsample
    L0_sub = datas$Lambda0_outer[subsample_index, subsample_index]
    R0_sub = cor_from_outer(L0_sub, param$sigma_sq_0)
    rm(datas); gc(verbose = FALSE)
    
    # compute BASIL
    fitBASIL = BASIL_point_estimates(Ys, Cs, k = param$k)
    # sample loadings given M = M_hat
    params_posterior_samples = compute_posterior_samples_cc(
      Ys, fitBASIL$Lambda_C, fitBASIL$Lambda_N, fitBASIL$tau_gamma,
      fitBASIL$tau_psi, fitBASIL$sigma_sq, fitBASIL$P_C)
    Lambda_outer_posterior_samples = sample_Lambda_outer(
      params_posterior_samples$Lambda_samples[subsample_index, , ])
    
    sigma_sq_hat = fitBASIL$sigma_sq
    rm(params_posterior_samples, fitBASIL, Ys); gc(verbose = FALSE)
    
    for (i in seq_len(dim(Lambda_outer_posterior_samples)[3])) {
      Lambda_outer_posterior_samples[, , i] =
        cor_from_outer(Lambda_outer_posterior_samples[, , i], sigma_sq_hat)
    }
    
    R_qs = apply(Lambda_outer_posterior_samples, c(1, 2),
                 function(x) quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
    rm(Lambda_outer_posterior_samples); gc(verbose = FALSE)
    
    # calculate coverage, off-diagonal only (the diagonal is 1 by construction)
    inside  = (R_qs[1, , ] < R0_sub) & (R_qs[2, , ] > R0_sub)
    offdiag = upper.tri(R0_sub) | lower.tri(R0_sub)
    
    data.frame(coverage = mean(inside[offdiag]),
               width = mean((R_qs[2, , ] - R_qs[1, , ])[offdiag]),
               zero_frac = zero_frac,
               sim = s,
               stringsAsFactors = FALSE)
  }
  
  res_list = mclapply(1:Nsim, one_sim, mc.cores = ncores, mc.preschedule = FALSE)
  
  coverage_df = do.call(rbind, res_list[ok])
  coverage_df$scenario = scenario_name
  coverage_df$p = as.character(param$p)
  rownames(coverage_df) = NULL
  
  cat(sprintf("%s: mean coverage %.4f (%d/%d simulations)\n",
              scenario_name, mean(coverage_df$coverage), sum(ok), Nsim))
  return(coverage_df)
}


#------------------------------------------------------------------------------#
#----------------------------# run simualations #------------------------------#


Nsim   = 25
ncores = 25   

# Setting 1: ~30% zero counts
param1_nb = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 60,
                 sd_gamma = 0.7, sd_psi = 0.1, mu_0 = 4, size_nb = 10)

# Setting 2: ~50% zero counts
param2_nb = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 60,
                 sd_gamma = 0.7, sd_psi = 0.1, mu_0 = -0.5, size_nb = 10)

#---# accuracy correlation #---#
df_sparse30_nb = run_simulation_studyNB(param1_nb, scenario_name = "sparsity30",
                                        Nsim = Nsim, seed = 463, ncores = ncores)

df_sparse50_nb = run_simulation_studyNB(param2_nb, scenario_name = "sparsity50",
                                        Nsim = Nsim, seed = 463, ncores = ncores)

Simboxplot_dfNB = rbind(df_sparse30_nb, df_sparse50_nb)
aggregate(cbind(err_cor, time, k_est, zero_frac) ~ model + scenario,
          data = Simboxplot_dfNB, FUN = mean)
aggregate(cbind(err_cor, time, k_est) ~ model + scenario,
          data = Simboxplot_dfNB, FUN = sd)


#---# uncertainty quantification #---#

subsample_size  = 200
subsample_index = 1:subsample_size
Nsimc  = 100

coverage_sparse30_nb = run_coverage_simulationNB(param1_nb, scenario_name = "sparsity30",
                                                 subsample_index, alpha = 0.05,
                                                 Nsim = Nsimc, seed = 463, ncores = ncores)

coverage_sparse50_nb = run_coverage_simulationNB(param2_nb, scenario_name = "sparsity50",
                                                 subsample_index, alpha = 0.05,
                                                 Nsim = Nsim, seed = 463, ncores = ncores)
SimUQ_nb = rbind(coverage_sparse30_nb, coverage_sparse50_nb)
aggregate(cbind(coverage, width, zero_frac) ~ scenario, data = SimUQ_nb, FUN = mean)
aggregate(cbind(coverage, width) ~ scenario, data = SimUQ_nb, FUN = sd)

