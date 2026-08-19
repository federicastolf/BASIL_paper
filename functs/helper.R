library(msigdbr)

source("helper/FACTOR_CODE_update.R")


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
  return(gene_set)
}


syntheticData = function(n, p, k, q, sigma_sq_0, sd_gamma, sd_psi, mseed, heteroscedastik=F){
  
  set.seed(mseed)
  C = get_geneSetMatrix(p, q)
  q = ncol(C)
  
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
  Psi_0 = Q_C %*% matrix(rnorm(p*k, 0, sd_psi), ncol=k)
  Lambda_0 = C_Gamma_0 + Psi_0
  Lambda0_outer = tcrossprod(Lambda_0)
  
  # data
  if(!heteroscedastik){
    Y = M_0 %*% t(Lambda_0) + sqrt(sigma_sq_0) * matrix(rnorm(n*p), nrow=n) 
  } else {
    Y = M_0 %*% t(Lambda_0) + matrix(rnorm(n*p), nrow=n) %*% diag(sqrt(runif(p, 0.5*sigma_sq_0, 2*sigma_sq_0)))
  }
  colnames(Y) = rownames(C)
  
  return(list("Y" = Y, "Lambda0_outer" = Lambda0_outer, "C" = C, 
              "M_0" = M_0))
}

syntheticDataPoisson = function(n, p, k, q, sd_gamma, sd_psi, mseed, offset = 2,
                                link = c("softplus", "exp"), nonneg_loadings = TRUE) {
  link = match.arg(link)
  set.seed(mseed)
  C = get_geneSetMatrix(p, q)
  q = ncol(C)
  
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
  M_0 = matrix(rnorm(n * k, 0, 1), ncol = k)
  
  # loadings
  Gamma_0 = matrix(rnorm(q * k, 0, sd_gamma), ncol = k)
  Psi_raw = matrix(rnorm(p * k, 0, sd_psi), ncol = k)
  if (nonneg_loadings) {
    Gamma_0 = abs(Gamma_0)
    Psi_raw = abs(Psi_raw)
  }
  C_Gamma_0 = C %*% Gamma_0
  Psi_0 = Q_C %*% Psi_raw   
  Lambda_0  = C_Gamma_0 + Psi_0
  Lambda0_outer = tcrossprod(Lambda_0)
  
  # link and rate  
  Z = M_0 %*% t(Lambda_0) 
  g = switch(link,
             softplus = function(x) log1p(exp(x)),
             exp = function(x) exp(x))
  mu = g(offset + Z)
  mu[mu < 1e-8] = 1e-8   # numerical safety
  
  # poisson
  Y = matrix(rpois(n * p, lambda = as.vector(mu)), nrow = n, ncol = p)
  colnames(Y) = rownames(C)
  
  Y_log = log1p(Y)
  Y_log = scale(Y_log, center = TRUE, scale = TRUE)
  attributes(Y_log)$`scaled:center` = NULL
  attributes(Y_log)$`scaled:scale` = NULL
  
  return(list(Y = Y, Y_log = Y_log, Lambda0_outer = Lambda0_outer, C = C, mu = mu))
}

# Negative binomial counts with the same factor structure on the loadings used in
# syntheticData(), imposed on the log-mean. 
syntheticDataNB = function(n, p, k, q, sigma_sq_0, sd_gamma, sd_psi, mseed,
                           mu_0 = 2, size_nb = 10) { 
  
  set.seed(mseed)
  C = get_geneSetMatrix(p, q)
  q = ncol(C)
  
  # compute P_C = C(C'C)^{-1}C' and Q_C = I_p - P_C
  tCCt = crossprod(C)
  s_cross_C = svd(tCCt)
  V_C = s_cross_C$u
  D_C_inv = diag(as.vector(1/sqrt(s_cross_C$d)))
  U_C = C %*% V_C %*% D_C_inv
  P_C = tcrossprod(U_C)
  Q_C = diag(1, p, p) - P_C
  
  # factors
  M_0 = matrix(rnorm(n*k, 0, 1), ncol = k)
  # loadings
  Gamma_0 = matrix(rnorm(q*k, 0, sd_gamma), ncol = k)
  C_Gamma_0 = C %*% Gamma_0
  Psi_0 = Q_C %*% matrix(rnorm(p*k, 0, sd_psi), ncol = k)
  Lambda_0 = C_Gamma_0 + Psi_0
  Lambda0_outer = tcrossprod(Lambda_0)
  
  # linear predictor on the log scale
  Eta = matrix(mu_0, nrow = n, ncol = p, byrow = TRUE) +
    M_0 %*% t(Lambda_0) +
    sqrt(sigma_sq_0) * matrix(rnorm(n*p), nrow = n)
  Mu = exp(pmin(Eta, eta_cap))
  
  # negative binomial counts
  Y_counts = matrix(rnbinom(n*p, size = size_nb, mu = as.vector(Mu)), nrow = n, ncol = p)
  Y_log = log1p(Y_counts)
  colnames(Y_counts) = colnames(Y_log) = colnames(Eta) = rownames(C)
  
  return(list("Y_counts" = Y_counts, "Y" = Y_log, "Eta" = Eta, "M_0" = M_0, 
              "Lambda_0" = Lambda_0, "Lambda0_outer" = Lambda0_outer, "C" = C))
}


run_simulation_study = function(param, scenario_name, Nsim, seed) {
  
  set.seed(seed)
  seeds_g = sample.int(9000, Nsim)
  
  # Initialize result storage
  results = list(kfitBASIL = numeric(Nsim), err_normBASIL = numeric(Nsim),
                 err_normBASIL_posterior = numeric(Nsim), 
                 err_factorsBASIL = numeric(Nsim), err_factorsBASIL_posterior = numeric(Nsim),
                 timeBASIL = numeric(Nsim), timeBASIL_posterior = numeric(Nsim), 
                 kfitROTATE = numeric(Nsim), err_normROTATE = numeric(Nsim), 
                 err_factorsROTATE = numeric(Nsim), timeROTATE = numeric(Nsim),
                 kfitPLIER = numeric(Nsim), err_normPLIER = numeric(Nsim), 
                 timePLIER = numeric(Nsim), 
                 err_normPLIER_ktrue = numeric(Nsim), timePLIER_ktrue = numeric(Nsim),
                 err_factorsPLIER_ktrue = numeric(Nsim),
                 parameters = param)
  
  if(is.null(param$heteroscedastik)){
    param$heteroscedastik = F
  }
  
  for (s in 1:Nsim) {
    # Simulate data
    datas = syntheticData(n = param$n, p = param$p, k = param$k, q = param$q,
                          sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                          sd_psi = param$sd_psi, mseed = seeds_g[s], heteroscedastik=param$heteroscedastik)
    Ys = datas$Y
    Cs = datas$C
    Lambda0_outer = datas$Lambda0_outer
    M_0 <- datas$M_0
    
    ## BASIL 
    ptmB <- proc.time()
    est_kBASIL <- estimate_latent_dimension(Ys, k_max = 50)
    fitBASIL <- BASIL_point_estimates(Ys, Cs, k = est_kBASIL$k_hat)
    etmB <- proc.time() - ptmB
    results$timeBASIL[s] <- etmB[1] + etmB[2]
    posterior_mean_BASIL <- compute_covariance_posterior_mean(Ys, fitBASIL)$Lambda_outer_mean
    etmB <- proc.time() - ptmB
    results$timeBASIL_posterior[s] <- etmB[1] + etmB[2]
    Lambda_BASIL <- fitBASIL$Lambda_C + fitBASIL$Lambda_N
    results$err_normBASIL[s] <- norm(tcrossprod(Lambda_BASIL) - Lambda0_outer, 
                                     type = 'F')/norm(Lambda0_outer, type = 'F')
    results$kfitBASIL[s] <- est_kBASIL$k_hat
    results$err_normBASIL_posterior[s] <- norm(posterior_mean_BASIL - Lambda0_outer,  
                                               type = 'F')/norm(Lambda0_outer, type = 'F')
    
    params_posterior_samples <- compute_posterior_samples_cc(
      Ys, fitBASIL$Lambda_C, fitBASIL$Lambda_N, fitBASIL$tau_gamma, 
      fitBASIL$tau_psi, fitBASIL$sigma_sq, fitBASIL$P_C)
    eta_samples <- sample_latent_factors(Ys, params_posterior_samples$Lambda_samples, params_posterior_samples$sigma_sq_samples)
    pr_hat <- procrustes(fitBASIL$M, M_0)
    results$err_factorsBASIL[s] <- norm(pr_hat$X.new - M_0, type = 'F')/norm(M_0, type = 'F')
    pr_tilde <- procrustes(eta_samples$M_mean, M_0)
    results$err_factorsBASIL_posterior[s] <-  norm(pr_tilde$X.new - M_0, type = 'F')/norm(M_0, type = 'F')
    
    
    
    ## ROTATE 
    K <- param$k
    startB <- matrix(rnorm(param$p * K), param$p, K)
    start <- list(B = startB, sigma = rep(1, K), theta = rep(0.5, K))
    ptmR <- proc.time()
    fitROTATE <- FACTOR_ROTATE(Y = Ys, lambda0 = 5, lambda1 = 0.001, start = start,
                               K = K, epsilon = 0.05, alpha = 1 / param$p, PX = TRUE, approximate = TRUE, 
                               stop = 100, varimax = TRUE, plot = FALSE)
    etmR <- proc.time() - ptmR
    results$timeROTATE[s] <- etmR[1] + etmR[2]
    Lambda_ROTATE <- fitROTATE$B
    results$err_normROTATE[s] <- norm(tcrossprod(Lambda_ROTATE) - Lambda0_outer,
                                      type = 'F')/norm(Lambda0_outer, type = 'F')
    M_hat_rotate <- latent_factor_full_conditional_mean(Ys, fitROTATE$B, fitROTATE$sigma^2)
    pr_rotate <- procrustes(M_hat_rotate, M_0)
    results$err_factorsROTATE[s]  <- norm(pr_rotate$X.new - M_0, type = 'F')/norm(M_0, type = 'F')
    
    ## PLIER 
    ptmP <- proc.time()
    fitPLIER <- PLIER(t(Ys), Cs, scale = FALSE, minGenes = 1, doCrossval = TRUE)
    etmP <- proc.time() - ptmP
    results$timePLIER[s] <- etmP[1] + etmP[2]
    results$kfitPLIER[s] <- nrow(fitPLIER$B)
    covPLIER <- (fitPLIER$Z %*% (fitPLIER$B %*% t(fitPLIER$B)) %*% t(fitPLIER$Z))/param$n
    results$err_normPLIER[s] <- norm(covPLIER - Lambda0_outer, type = 'F')/norm(Lambda0_outer, type = 'F')
    ## PLIER (true k)
    ptmP <- proc.time()
    fitPLIER <- PLIER(t(Ys), Cs, scale = FALSE, minGenes = 1, doCrossval = TRUE, k=param$k)
    etmP <- proc.time() - ptmP
    results$timePLIER_ktrue[s] <- etmP[1] + etmP[2]
    covPLIER <- (fitPLIER$Z %*% (fitPLIER$B %*% t(fitPLIER$B)) %*% t(fitPLIER$Z))/param$n
    results$err_normPLIER_ktrue[s] <- norm(covPLIER - Lambda0_outer, type = 'F')/norm(Lambda0_outer, type = 'F')
    s_plier <- svd(t(fitPLIER$B))
    pr_plier <- procrustes(sqrt(param$n) * s_plier$u, M_0)
    results$err_factorsPLIER_ktrue[s] <- norm(pr_plier$X.new - M_0, type = 'F')/norm(M_0, type = 'F')
    
    
    if (s %% 5 == 0) {
      cat(sprintf("Completed simulation %d/%d\n", s, Nsim))
    }
  }
  
  #Convert to tidy data frame
  df_basil <- data.frame(err_norm = results$err_normBASIL, time = results$timeBASIL,
                         k_est = results$kfitBASIL, err_factors = results$err_factorsBASIL,
                         model = "BASIL", p = param$p, scenario = scenario_name,
                         stringsAsFactors = FALSE)
  
  df_basil_post <- data.frame(err_norm = results$err_normBASIL_posterior,
                              time = results$timeBASIL_posterior, k_est = results$kfitBASIL,
                              err_factors = results$err_factorsBASIL_posterior, model = "BASIL_posterior",
                              p = param$p,scenario = scenario_name, stringsAsFactors = FALSE)
  
  df_rotate <- data.frame(err_norm = results$err_normROTATE, time = results$timeROTATE,
                          err_factors = results$err_factorsROTATE,
                          k_est = results$kfitROTATE, model = "ROTATE", p = param$p, scenario = scenario_name, 
                          stringsAsFactors = FALSE)
  
  df_plier <- data.frame(err_norm = results$err_normPLIER, time = results$timePLIER,
                         k_est = results$kfitPLIER,
                         model = "PLIER", p = param$p, scenario = scenario_name, 
                         err_factors = rep(NA, Nsim),
                         stringsAsFactors = FALSE)
  
  df_plier_ktrue <- data.frame(err_norm = results$err_normPLIER_ktrue, time = results$timePLIER_ktrue, err_factors = results$err_factorsPLIER_ktrue,
                               k_est = param$k, model = "PLIER-true k", p = param$p, scenario = scenario_name,
                               stringsAsFactors = FALSE)
  
  # Combine all methods
  df_combined <- rbind(df_basil, df_basil_post, df_rotate, df_plier, df_plier_ktrue)
  
  return(df_combined)
}
# Function to run coverage simulation
run_coverage_simulation <- function(param, scenario_name, subsample_index, alpha, 
                                    Nsim, seed){
  
  set.seed(seed)
  seeds_g <- sample.int(9000, Nsim)
  ccCoverage <- numeric(Nsim)
  
  # Run simulations
  for (s in 1:Nsim) {
    # Simulate data
    datas <- syntheticData(n = param$n, p = param$p, k = param$k, q = param$q,
                           sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                           sd_psi = param$sd_psi, mseed = seeds_g[s])
    Ys <- datas$Y
    Cs <- datas$C
    Lambda0_outer <- datas$Lambda0_outer
    
    # Compute BASIL
    fitBASIL <- BASIL_point_estimates(Ys, Cs, k = param$k)
    # Sample loadings given M = M_hat
    params_posterior_samples <- compute_posterior_samples_cc(
      Ys, fitBASIL$Lambda_C, fitBASIL$Lambda_N, fitBASIL$tau_gamma, 
      fitBASIL$tau_psi, fitBASIL$sigma_sq, fitBASIL$P_C)
    Lambda_outer_posterior_samples <- sample_Lambda_outer(
      params_posterior_samples$Lambda_samples[subsample_index, , ])
    Lambda_outer_qs <- apply(Lambda_outer_posterior_samples, c(1, 2),
                             function(x) quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
    
    # Calculate coverage
    cov <- mean((Lambda_outer_qs[1, , ] < Lambda0_outer[subsample_index, subsample_index]) &
                  (Lambda_outer_qs[2, , ] > Lambda0_outer[subsample_index, subsample_index]))
    ccCoverage[s] <- cov
    
    if (s %% 5 == 0) {
      cat(sprintf("Completed simulation %d/%d, current mean coverage: %.4f\n", 
                  s, Nsim, mean(ccCoverage[1:s])))
    }
  }
  coverage_df = data.frame(coverage = ccCoverage, scenario = scenario_name, 
                           p = as.character(param$p), stringsAsFactors = FALSE)
  
  return(coverage_df)
}


# Function to create scatter plot with identity line
plot_correlation_scatter = function(data, title, lim_ax, point_color = "#1170aa") {
  ggplot(data, aes(x = Observed, y = Predicted)) +
    geom_point(color = point_color, shape = 1) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1) +
    labs(title = title, x = "Observed", y = "Predicted") +
    xlim(lim_ax) + ylim(lim_ax) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 16))
}

latent_factor_full_conditional_mean <- function(Y, Lambda, Sigma){
  k <- ncol(Lambda)
  mean_f <- Y %*% diag(1/Sigma) %*% Lambda %*% solve(t(Lambda) %*% diag(1/Sigma) %*% Lambda + diag(1, k, k))
  return(mean_f)
}
