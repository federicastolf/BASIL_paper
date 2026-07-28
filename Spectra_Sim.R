library(BASIL)
library(tidyverse)

rm(list=ls())

source("helper.R")
Nsim = 25

#--------------# simulate data and save for export in Python #-----------------#

set.seed(123)
seed_v = sample.int(900, Nsim)

data_dir = "sim_data"
dir.create(data_dir, showWarnings = FALSE)
data_poin = vector("list", Nsim)

for (s in 1:Nsim) {
  # positive loadings
  data_poin[[s]] = syntheticDataPoisson(n = 500, p = 1000, k = 10, q = 500,
                                        sd_gamma = 0.7, sd_psi = 0.1, 
                                        mseed = seed_v[s], offset = 2, 
                                        link = "softplus", nonneg_loadings = T)
  
  # mixed-sign loadings
  # data_poin[[s]] = syntheticDataPoisson(n = 500, p = 1000, k = 10, q = 500,
  #                                       sd_gamma = 0.7, sd_psi = 0.1, 
  #                                       mseed = seed_v[s], offset = 2, 
  #                                       link = "softplus", nonneg_loadings = F)
  
  saveRDS(list(Y = data_poin[[s]]$Y, C = data_poin[[s]]$C,
               Lambda0_outer = data_poin[[s]]$Lambda0_outer),
          file = file.path(data_dir, sprintf("sim_%03d.rds", s)))
}

for (s in 1:Nsim) {
  d = file.path(data_dir, sprintf("sim_%03d", s))
  dir.create(d, showWarnings = FALSE)
  write.table(data_poin[[s]]$Y,             file.path(d, "Y.csv"),
              sep=",", row.names=FALSE, col.names=FALSE)
  write.table(data_poin[[s]]$C,             file.path(d, "C.csv"),
              sep=",", row.names=FALSE, col.names=FALSE)
  write.table(data_poin[[s]]$Lambda0_outer, file.path(d, "Lambda0_outer.csv"),
              sep=",", row.names=FALSE, col.names=FALSE)
  writeLines(colnames(data_poin[[s]]$C), file.path(d, "gene_set_names.txt"))
}

#---------# load Spectra results from Python and compare with BASIL #----------#

results = list(kfitBASIL = numeric(Nsim), err_normBASIL = numeric(Nsim),
               timeBASIL = numeric(Nsim), kfitSpectra = numeric(Nsim),
               err_normSpectra = numeric(Nsim), timeSpectra = numeric(Nsim))
spectra_dir = "spectra_results"

for (s in 1:Nsim) {
  
  datas = data_poin[[s]]
  Ys = datas$Y_log
  Cs = datas$C
  Lambda0_outer = datas$Lambda0_outer
  n = nrow(Ys)
  
  ## ---- BASIL ----
  est_kBASIL = estimate_latent_dimension(Ys, k_max = 50)
  ptmB = Sys.time()
  fitBASIL   = BASIL_point_estimates(Ys, Cs, k = est_kBASIL$k_hat)
  etmB = Sys.time() - ptmB
  results$timeBASIL[s] = as.numeric(etmB, units = "secs")
  results$kfitBASIL[s] = est_kBASIL$k_hat
  Lambda_BASIL = fitBASIL$Lambda_C + fitBASIL$Lambda_N
  cov_BASIL = tcrossprod(Lambda_BASIL)
  results$err_normBASIL[s] = norm(cov_BASIL - Lambda0_outer, "F") /
    norm(Lambda0_outer, "F")

  ## ---- load results from SPECTRA ----
  sim_id = sprintf("sim_%03d", s)
  sim_folder = file.path(spectra_dir, sim_id)

  Lambda_S = as.matrix(read.csv(file.path(sim_folder, "Lambda_hat.csv"), 
                                header = FALSE))
  M_S = as.matrix(read.csv(file.path(sim_folder, "M_hat.csv"), header = FALSE))
  meta = readLines(file.path(sim_folder, "meta.txt"))
  time_S = as.numeric(strsplit(meta[grep("^time_seconds", meta)], "\t")[[1]][2])
  K_S = as.integer(strsplit(meta[grep("^K", meta)], "\t")[[1]][2])
  
  # low-rank covariance for SPECTRA: (Lambda M'M Lambda') / n
  MtM_over_n = crossprod(M_S) / n
  cov_Spectra = Lambda_S %*% MtM_over_n %*% t(Lambda_S)
  results$timeSpectra[s] = time_S
  results$kfitSpectra[s] = K_S
  results$err_normSpectra[s] = norm(cov_Spectra - Lambda0_outer, "F") /
    norm(Lambda0_outer, "F")
  
}

