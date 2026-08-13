library(LaplacesDemon)
library(BASIL)
rm(list = ls())
source("functs/face.R")
source("functs/helper.R")

# high
param_h = list(n = 100, p = 300, k = 10, q = 150, sigma_sq_0 = 15, sd_gamma = 0.7,
               sd_psi = 0.1)
# low
param_l = list(n = 100, p = 300, k = 10, q = 150, sigma_sq_0 = 15, sd_gamma = 0.4,
               sd_psi = 0.7)

Nsim = 25
base_seed = 433
subsample_index = 1:100
alpha_cov = 0.05

run_simFACE = function(s, param, subsample_index, alpha_cov) {
  p = param$p
  datas = syntheticData(n = param$n, p = param$p, k = param$k, q = param$q,
                        sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                        sd_psi = param$sd_psi, mseed = base_seed + s)
  Ys = datas$Y
  Cs = datas$C
  Lambda0_outer = datas$Lambda0_outer
  true_sub = Lambda0_outer[subsample_index, subsample_index]
  
  ## ---- BASIL ----
  ptmB = proc.time()
  est_kBASIL = estimate_latent_dimension(Ys, k_max = 50)
  fitBASIL = BASIL_point_estimates(Ys, Cs, k = est_kBASIL$k_hat)
  params_posterior_samples = compute_posterior_samples_cc(
    Ys, fitBASIL$Lambda_C, fitBASIL$Lambda_N, fitBASIL$tau_gamma,
    fitBASIL$tau_psi, fitBASIL$sigma_sq, fitBASIL$P_C)
  etmB = proc.time() - ptmB
  Lambda_BASIL = fitBASIL$Lambda_C + fitBASIL$Lambda_N
  err_BASIL = norm(tcrossprod(Lambda_BASIL) - Lambda0_outer, type = "F") /
    norm(Lambda0_outer, type = "F")
  # coverage
  Lambda_outer_posterior_samples = sample_Lambda_outer(
    params_posterior_samples$Lambda_samples[subsample_index, , ])
  Lambda_outer_qs_BASIL = apply(
    Lambda_outer_posterior_samples, c(1, 2),
    function(x) quantile(x, probs = c(alpha_cov / 2, 1 - alpha_cov / 2)))
  cov_BASIL = mean((Lambda_outer_qs_BASIL[1, , ] < true_sub) &
                     (Lambda_outer_qs_BASIL[2, , ] > true_sub))
  
  ## ---- FACE (CMR_cusp_GS) ----
  kF = floor(p / 2)
  ptmF = proc.time()
  out.cmr.cusp = CMR_cusp_GS(Ys, Cs, k = kF, S = 10000, burnin = 5000,
                             my.seed = 100 + s, alpha = 5, a.theta = 1/2,
                             b.theta = 1/2, which.cov.group = NA)
  etmF = proc.time() - ptmF
  Lambda.mean = matrix(colMeans(out.cmr.cusp$Lambda), nrow = p, ncol = kF)
  err_FACE = norm(tcrossprod(Lambda.mean) - Lambda0_outer, type = "F") /
    norm(Lambda0_outer, type = "F")
  
  # coverage
  n.draws = nrow(out.cmr.cusp$cov)
  sub.draws = sapply(1:n.draws, function(i) {
    cov.mat = matrix(out.cmr.cusp$cov[i, ], nrow = p, ncol = p)
    cov.mat[subsample_index, subsample_index]
  }, simplify = "array")
  Lambda_outer_qs_FACE = apply(
    sub.draws, c(1, 2),
    function(x) quantile(x, probs = c(alpha_cov / 2, 1 - alpha_cov / 2)))
  cov_FACE = mean((Lambda_outer_qs_FACE[1, , ] < true_sub) &
                    (Lambda_outer_qs_FACE[2, , ] > true_sub))
  
  list(replica = s, err_normBASIL = err_BASIL, timeBASIL = unname(etmB["elapsed"]),
       covBASIL = cov_BASIL, err_normFACE = err_FACE,
       timeFACE = unname(etmF["elapsed"]), covFACE = cov_FACE)
}

res_list_h = vector("list", Nsim)
res_list_l = vector("list", Nsim)

for (s in 1:Nsim) {
  res_list_h[[s]] = run_simFACE(s, param_h, subsample_index, alpha_cov)
}

for (s in 1:Nsim) {
  res_list_l[[s]] = run_simFACE(s, param_l, subsample_index, alpha_cov)
}