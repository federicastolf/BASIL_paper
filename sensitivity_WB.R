library(BASIL)
library(PLIER)           
library(tidyverse)
library(reshape2)
library(gridExtra)
#library(ComplexHeatmap)
library(scico)
#library(forcats)
library(tidytext)

rm(list=ls())

#----------# Load whole blood gene expression data #-----------#

data("dataWholeBlood", package = "PLIER") 
data("bloodCellMarkersIRISDMAP", package = "PLIER")
data("canonicalPathways", package = "PLIER")
allPaths = PLIER::combinePaths(bloodCellMarkersIRISDMAP, canonicalPathways)
cm = intersect(rownames(dataWholeBlood), rownames(allPaths))
allPaths = allPaths[cm,]
dataWholeBlood = dataWholeBlood[cm,]

Y_full = t(dataWholeBlood)
Y_full = scale(Y_full)
n = nrow(Y_full)
p = ncol(Y_full)

# select train test
n_train = round(0.8 * n)
n_test  = n - n_train   # ~20% of the samples held out

# Two different random training subsets 
set.seed(213)
idx_train_cfg1 = sample(seq_len(n), n_train)
set.seed(132)
idx_train_cfg2 = sample(seq_len(n), n_train)
train_idx_list = list("cfg1" = idx_train_cfg1, "cfg2" = idx_train_cfg2)

#------# Fit BASIL on each training set #---------#

fits = list()
for (cfg in names(train_idx_list)) {
  Y_train = Y_full[train_idx_list[[cfg]], ]
  est_k = estimate_latent_dimension(Y_train, k_max = 12)
  fits[[cfg]] = BASIL_point_estimates(Y = Y_train, C = allPaths, k = est_k$k_hat,
                                       v0 = 1, sigma_sq0 = 1)
}
sapply(fits, function(f) f$tau_gamma / f$tau_psi)

# Top 100 genes by empirical variance
gene_variances = apply(t(dataWholeBlood), 2, var)
idx_genes = order(gene_variances, decreasing = TRUE)[1:100]

#-------# compute correlation + UQ #----------#
compute_basil_corr_zeroed = function(fit, Y, idx, v0 = 1, sigma_sq_0 = 1,
                                      n_MC = 500, ci_level = 0.95) {
  
  loadings_samples = compute_posterior_samples_cc(
    Y = Y, Lambda_C = fit$Lambda_C, Lambda_N = fit$Lambda_N, tau_gamma = fit$tau_gamma,
    tau_psi = fit$tau_psi, sigma_sq = fit$sigma_sq, P_C = fit$P_C, v0 = v0, 
    sigma_sq_0 = sigma_sq_0, n_MC = n_MC
  )
  cor_posterior = compute_correlation_posterior_samples_cc(
    Lambda_samples = loadings_samples$Lambda_samples[idx, TRUE, TRUE],
    sigma_sq_samples = loadings_samples$sigma_sq_samples, samples = TRUE
  )
  
  corr_mean = cor_posterior$posterior_mean
  corr_mean_zeroed = corr_mean
  alpha = 1 - ci_level
  corr_qs = apply(cor_posterior$posterior_samples, c(1, 2), quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))
  corr_mean_zeroed[(corr_qs[1, , ] < 0) & (corr_qs[2, , ] > 0)] <- 0
  
  list(loadings_samples  = loadings_samples, cor_posterior = cor_posterior,
       corr_mean = corr_mean, corr_mean_zeroed = corr_mean_zeroed,
       prop_zeroed = mean(corr_mean_zeroed == 0))
}

results = Map(function(fit, cfg) {
  Y_train = Y_full[train_idx_list[[cfg]], ]
  compute_basil_corr_zeroed(fit, Y = Y_train, idx = idx_genes)
}, fits, names(fits))

gene_names = colnames(Y_full)[idx_genes]
for (cfg in names(results)) {
  rownames(results[[cfg]]$corr_mean_zeroed) = gene_names
  colnames(results[[cfg]]$corr_mean_zeroed) = gene_names
  rownames(results[[cfg]]$corr_mean) = gene_names
  colnames(results[[cfg]]$corr_mean) = gene_names
}
print(sapply(results, function(x) x$prop_zeroed))

# setting for plot
cols = scico(51, palette = "vik") 
breaks = seq(-1, 1, length.out = 52)
Corr_empirical = cor(Y_full[, idx_genes])
emp_corr_plot = plot_correlation_heatmap(Corr_empirical, title = "Empirical", 
                                          cluster = TRUE, breaks = breaks)
row_order = emp_corr_plot$tree_row$order
col_order = emp_corr_plot$tree_col$order

basilUQ_corr_plot_cfg1 = plot_correlation_heatmap(
  results[["cfg1"]]$corr_mean_zeroed, title = "", row_order = row_order, 
  col_order = col_order, breaks = breaks
)
basilUQ_corr_plot_cfg2 = plot_correlation_heatmap(
  results[["cfg2"]]$corr_mean_zeroed, title = "", row_order = row_order, 
  col_order = col_order, breaks = breaks
)
basil_corr_plot_cfg1 = plot_correlation_heatmap(
  results[["cfg1"]]$corr_mean, title = "", row_order = row_order, 
  col_order = col_order, breaks = breaks
)
basil_corr_plot_cfg2 <- plot_correlation_heatmap(
  results[["cfg2"]]$corr_mean, title = "", row_order = row_order, 
  col_order = col_order, breaks = breaks
)

corrWB_cfg = grid.arrange(
  basilUQ_corr_plot_cfg1$gtable, basilUQ_corr_plot_cfg2$gtable, 
  basil_corr_plot_cfg1$gtable, basil_corr_plot_cfg2$gtable, ncol = 2
)

# ggsave('CorrWB_sens_traintest.jpeg', plot = corrWB_cfg, width = 6.5, height = 5.8)


#--------# posterior samples of Gamma #--------#

# top pathways for the first three factors
p_names = colnames(allPaths)
idf_fa = 1:3

gamma_results = list()
for (cfg in names(results)) {
  Gamma_postsamples = compute_Gamma_samples(
    results[[cfg]]$loadings_samples$Lambda_samples, allPaths
  )
  ci_lowerG = apply(Gamma_postsamples, c(1, 2), quantile, probs = 0.025)
  ci_upperG = apply(Gamma_postsamples, c(1, 2), quantile, probs = 0.975)
  contains_zeroG = (ci_lowerG < 0) & (ci_upperG > 0)
  print(table(contains_zeroG) / length(contains_zeroG))
  Gamma_UQ = apply(Gamma_postsamples, c(1, 2), mean)
  Gamma_UQ[contains_zeroG] = 0
  
  data_pplot_cfg = get_topPathways(Gamma_UQ[, idf_fa], p_names, top_n = 5,
                                    clean_names = TRUE, ci_lower = ci_lowerG[, idf_fa],
                                    ci_upper = ci_upperG[, idf_fa], max_length = 90)
  data_pplot_cfg = data_pplot_cfg[, 1:3]
  data_pplot_cfg$config = cfg
  gamma_results[[cfg]] = data_pplot_cfg
}

data_combined = bind_rows(gamma_results) %>%
  mutate(AbsLoading = abs(Loading), Pathway_wrapped = str_wrap(Pathway, width = 30))

make_factor_plot = function(f_num) {
  df = data_combined %>% filter(Factor_num == f_num)
  
  ord = df %>%
    group_by(Pathway_wrapped) %>%
    summarise(order_val = max(AbsLoading), .groups = "drop") %>%
    arrange(desc(order_val)) %>%
    pull(Pathway_wrapped)
  
  df = df %>%
    mutate(rank = match(Pathway_wrapped, ord),
           y_num = rank + ifelse(config == "cfg2", 0.12, -0.12),
           pt_stroke = ifelse(config == "cfg2", 1.6, 1))
  
  ggplot(df, aes(x = AbsLoading, y = y_num, shape = config, stroke = pt_stroke)) +
    geom_point(size = 2.5, color = "black") +
    scale_shape_manual(values = c("cfg1" = 16, "cfg2" = 4)) +
    scale_y_continuous(breaks = seq_along(ord), labels = ord, trans = "reverse") +
    xlim(c(0, 1.6)) +
    labs(x = "Absolute loading", y = NULL, title = paste0("Factor ", f_num)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.position = "none",
          axis.text.y = element_text(size = 8, lineheight = 0.8))
}

plots = lapply(1:3, make_factor_plot)
fplot_sens_traintest = grid.arrange(grobs = plots, nrow = 1)

# ggsave('sens_dotplot_train.jpeg', plot = fplot_sens_traintest, width = 11, height = 4)



###################################################################################
###################################################################################

rm(list=ls())

# Load whole blood gene expression data
data("dataWholeBlood", package = "PLIER") 
data("bloodCellMarkersIRISDMAP", package = "PLIER")
data("canonicalPathways", package = "PLIER")
allPaths = PLIER::combinePaths(bloodCellMarkersIRISDMAP, canonicalPathways)
cm = intersect(rownames(dataWholeBlood), rownames(allPaths))
allPaths = allPaths[cm,]
dataWholeBlood = dataWholeBlood[cm,]
dim(dataWholeBlood)  
dim(allPaths)        

Y = t(dataWholeBlood)
Y = scale(Y)
n = nrow(Y)
p = ncol(Y)

#------# fit BASIL for different k #--------#

fit_basil8 = BASIL_point_estimates(Y = Y, C = allPaths, k = 8, v0 = 1, sigma_sq0 = 1)
fit_basil5 = BASIL_point_estimates(Y = Y, C = allPaths, k = 5, v0 = 1, sigma_sq0 = 1)
fit_basil12 = BASIL_point_estimates(Y = Y, C = allPaths, k = 12, v0 = 1, sigma_sq0 = 1)

fit_basil8$tau_gamma / fit_basil8$tau_psi
fit_basil5$tau_gamma / fit_basil5$tau_psi
fit_basil12$tau_gamma / fit_basil12$tau_psi

# Top 100 genes by empirical variance
gene_variances = apply(t(dataWholeBlood), 2, var)
idx = order(gene_variances, decreasing = TRUE)[1:100]
Y_subset = Y[,idx]

#-------# compute correlation + UQ #----------#

compute_basil_corr_zeroed = function(fit, Y, idx, v0 = 1, sigma_sq_0 = 1,
                                      n_MC = 500, ci_level = 0.95) {
  
  loadings_samples = compute_posterior_samples_cc(
    Y = Y, Lambda_C = fit$Lambda_C, Lambda_N = fit$Lambda_N, tau_gamma = fit$tau_gamma,
    tau_psi = fit$tau_psi, sigma_sq = fit$sigma_sq, P_C = fit$P_C, v0 = v0, 
    sigma_sq_0 = sigma_sq_0, n_MC = n_MC
  )
  cor_posterior = compute_correlation_posterior_samples_cc(
    Lambda_samples = loadings_samples$Lambda_samples[idx, TRUE, TRUE],
    sigma_sq_samples = loadings_samples$sigma_sq_samples, samples = TRUE
  )
  
  corr_mean = cor_posterior$posterior_mean
  corr_mean_zeroed = corr_mean
  alpha = 1 - ci_level
  corr_qs = apply(cor_posterior$posterior_samples, c(1, 2), quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))
  corr_mean_zeroed[(corr_qs[1, , ] < 0) & (corr_qs[2, , ] > 0)] <- 0
  
  list(loadings_samples  = loadings_samples, cor_posterior = cor_posterior,
       corr_mean = corr_mean, corr_mean_zeroed = corr_mean_zeroed,
       prop_zeroed = mean(corr_mean_zeroed == 0))
}
fits = list("5" = fit_basil5, "8" = fit_basil8, "12" = fit_basil12)
results = lapply(fits, compute_basil_corr_zeroed, Y = Y, idx = idx)

for (nm in names(results)) {
  rownames(results[[nm]]$corr_mean_zeroed) = colnames(Y_subset)
  colnames(results[[nm]]$corr_mean_zeroed) = colnames(Y_subset)
  rownames(results[[nm]]$corr_mean) = colnames(Y_subset)
  colnames(results[[nm]]$corr_mean) = colnames(Y_subset)
}
print(sapply(results, function(x) x$prop_zeroed))

# setting for plot
cols = scico(51, palette = "vik") 
breaks = seq(-1, 1, length.out = 52)

Corr_empirical = cor(Y_subset)
emp_corr_plot = plot_correlation_heatmap(Corr_empirical, title = "Empirical", 
                                          cluster = TRUE, breaks = breaks)
row_order = emp_corr_plot$tree_row$order
col_order = emp_corr_plot$tree_col$order

basilUQ_corr_plot8 = plot_correlation_heatmap(
  results[["8"]]$corr_mean_zeroed, title = "k=8", row_order = row_order, 
  col_order = col_order, breaks = breaks
)
basilUQ_corr_plot5 = plot_correlation_heatmap(
  results[["5"]]$corr_mean_zeroed, title = "k=5", row_order = row_order, 
  col_order = col_order, breaks = breaks
)
basilUQ_corr_plot12 = plot_correlation_heatmap(
  results[["12"]]$corr_mean_zeroed, title = "k=12", row_order = row_order, 
  col_order = col_order, breaks = breaks
)
basil_corr_plot8 =  plot_correlation_heatmap(
  results[["8"]]$corr_mean, title = "k=8", row_order = row_order, 
  col_order = col_order, breaks = breaks
)
basil_corr_plot5 = plot_correlation_heatmap(
  results[["5"]]$corr_mean, title = "k=5", row_order = row_order, 
  col_order = col_order, breaks = breaks
)
basil_corr_plot12  = plot_correlation_heatmap(
  results[["12"]]$corr_mean, title = "k=12", row_order = row_order, 
  col_order = col_order, breaks = breaks
)

corrWB = grid.arrange(basilUQ_corr_plot5$gtable, basilUQ_corr_plot8$gtable, 
                       basilUQ_corr_plot12$gtable, basil_corr_plot5$gtable,
                       basil_corr_plot8$gtable, basil_corr_plot12$gtable, ncol = 3)

# ggsave('CorrWB_sens_k.jpeg', plot=corrWB, width = 9.1, height = 5.8)


#-----# posterior samples of Gamma #-------#
p_names = colnames(allPaths)
idf_fa = c(1:3)

# 5
Gamma_postsamples = compute_Gamma_samples(results[["5"]]$loadings_samples$Lambda_samples,
                                           allPaths)
ci_lowerG = apply(Gamma_postsamples, c(1, 2), quantile, probs = 0.025)
ci_upperG = apply(Gamma_postsamples, c(1, 2), quantile, probs = 0.975)
contains_zeroG = (ci_lowerG < 0) & (ci_upperG > 0)
table(contains_zeroG)/length(contains_zeroG) 
Gamma_UQ = apply(Gamma_postsamples, c(1, 2), mean)
Gamma_UQ[contains_zeroG] = 0
data_pplot = get_topPathways(Gamma_UQ[,idf_fa], p_names, top_n = 5,
                              clean_names = TRUE, ci_lower = ci_lowerG[,idf_fa],
                              ci_upper = ci_upperG[,idf_fa], max_length = 90)

# 8
Gamma_postsamples8 = compute_Gamma_samples(results[["8"]]$loadings_samples$Lambda_samples,
                                           allPaths)
ci_lowerG8 = apply(Gamma_postsamples8, c(1, 2), quantile, probs = 0.025)
ci_upperG8 = apply(Gamma_postsamples8, c(1, 2), quantile, probs = 0.975)
contains_zeroG8 = (ci_lowerG8 < 0) & (ci_upperG8 > 0)
table(contains_zeroG8)/length(contains_zeroG8) 
Gamma_UQ8 = apply(Gamma_postsamples8, c(1, 2), mean)
Gamma_UQ8[contains_zeroG8] = 0
data_pplot8 = get_topPathways(Gamma_UQ8[,idf_fa], p_names, top_n = 5,
                              clean_names = TRUE, ci_lower = ci_lowerG8[,idf_fa],
                              ci_upper = ci_upperG8[,idf_fa], max_length = 90)

# 12
Gamma_postsamples12 = compute_Gamma_samples(results[["12"]]$loadings_samples$Lambda_samples,
                                            allPaths)
ci_lowerG12 = apply(Gamma_postsamples12, c(1, 2), quantile, probs = 0.025)
ci_upperG12 = apply(Gamma_postsamples12, c(1, 2), quantile, probs = 0.975)
contains_zeroG12 = (ci_lowerG12 < 0) & (ci_upperG12 > 0)
table(contains_zeroG12)/length(contains_zeroG12) 
Gamma_UQ12 = apply(Gamma_postsamples12, c(1, 2), mean)
Gamma_UQ12[contains_zeroG12] = 0
data_pplot12 = get_topPathways(Gamma_UQ12[,idf_fa], p_names, top_n = 5,
                               clean_names = TRUE, ci_lower = ci_lowerG12[,idf_fa],
                               ci_upper = ci_upperG12[,idf_fa], max_length = 90)

# data for plot
data_pplot12 = data_pplot12[,1:3]
data_pplot8 = data_pplot8[,1:3]
data_pplot5 = data_pplot[,1:3]

data_pplot5$k = 5
data_pplot8$k = 8
data_pplot12$k =  12

data_combined = bind_rows(data_pplot5, data_pplot8, data_pplot12)

data_combined = data_combined %>%
  mutate(AbsLoading = abs(Loading),
         Pathway_wrapped = str_wrap(Pathway, width = 30))

make_factor_plot = function(f_num) {
  df = data_combined %>% filter(Factor_num == f_num)
  
  ord = df %>%
    group_by(Pathway_wrapped) %>%
    summarise(order_val = max(AbsLoading), .groups = "drop") %>%
    arrange(desc(order_val)) %>%
    pull(Pathway_wrapped)
  
  df = df %>%
    mutate(rank = match(Pathway_wrapped, ord),
           y_num = rank + case_when(k == 5  ~ -0.2,
                                    k == 8  ~   0,
                                    k == 12 ~  0.2),
           pt_stroke = ifelse(k == 5, 1.2, 0.8))  
  
  ggplot(df, aes(x = AbsLoading, y = y_num, shape = factor(k), stroke = pt_stroke)) +
    geom_point(size = 2.5, color = "black") +
    scale_shape_manual(values = c("5" = 4, "8" = 16, "12" = 15)) +
    scale_y_continuous(breaks = seq_along(ord), labels = ord, trans = "reverse") +
    xlim(0, 1.5) +
    labs(x = "Absolute loading", y = NULL, shape = "k",
         title = paste0("Factor ", f_num)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.position = "none",
          axis.text.y = element_text(size = 8, lineheight = 0.8))
}

plots = lapply(1:3, make_factor_plot)
fplot_sens = grid.arrange(grobs = plots, nrow = 1)

# ggsave('sens_dotplot_k.jpeg', plot = fplot_sens, width = 11, height = 4)




###################################################################################
###################################################################################

rm(list = ls())

data("dataWholeBlood",             package = "PLIER")
data("bloodCellMarkersIRISDMAP",   package = "PLIER")
data("canonicalPathways",          package = "PLIER")
allPaths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, canonicalPathways)
cm       <- intersect(rownames(dataWholeBlood), rownames(allPaths))
allPaths <- allPaths[cm, ]
dataWholeBlood <- dataWholeBlood[cm, ]

Y <- scale(t(dataWholeBlood))    
n <- nrow(Y)
p <- ncol(Y)

# Blom scores: qnorm((r - 3/8) / (n + 1/4)).
int_blom <- function(x) {
  ok <- !is.na(x)
  z  <- rep(NA_real_, length(x))
  r  <- rank(x[ok], ties.method = "average")
  z[ok] <- qnorm((r - 3/8) / (sum(ok) + 1/4))
  z
}

Y1 <- apply(Y, 2, int_blom)                 

fit_orig <- BASIL_point_estimates(Y = Y,  C = allPaths, k = 8, v0 = 1, sigma_sq0 = 1)
fit_int  <- BASIL_point_estimates(Y = Y1, C = allPaths, k = 8, v0 = 1, sigma_sq0 = 1)

fits <- list(Original = fit_orig, `Rank-INT` = fit_int)
Ys   <- list(Original = Y,        `Rank-INT` = Y1)

c(Original = fit_orig$tau_gamma / fit_orig$tau_psi,
  `Rank-INT` = fit_int$tau_gamma / fit_int$tau_psi)

L_orig <- fit_orig$Lambda_C + fit_orig$Lambda_N
L_int  <- fit_int$Lambda_C  + fit_int$Lambda_N

gene_variances <- apply(t(dataWholeBlood), 2, var)
idx      <- order(gene_variances, decreasing = TRUE)[1:100]
Y_subset <- Y[, idx]

compute_basil_corr_zeroed <- function(fit, Y, idx, v0 = 1, sigma_sq_0 = 1,
                                      n_MC = 500, ci_level = 0.95) {
  
  loadings_samples <- compute_posterior_samples_cc(
    Y = Y, Lambda_C = fit$Lambda_C, Lambda_N = fit$Lambda_N,
    tau_gamma = fit$tau_gamma, tau_psi = fit$tau_psi, sigma_sq = fit$sigma_sq,
    P_C = fit$P_C, v0 = v0, sigma_sq_0 = sigma_sq_0, n_MC = n_MC
  )
  cor_posterior <- compute_correlation_posterior_samples_cc(
    Lambda_samples   = loadings_samples$Lambda_samples[idx, TRUE, TRUE],
    sigma_sq_samples = loadings_samples$sigma_sq_samples, samples = TRUE
  )
  
  corr_mean        <- cor_posterior$posterior_mean
  corr_mean_zeroed <- corr_mean
  alpha   <- 1 - ci_level
  corr_qs <- apply(cor_posterior$posterior_samples, c(1, 2), quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))
  corr_mean_zeroed[(corr_qs[1, , ] < 0) & (corr_qs[2, , ] > 0)] <- 0
  
  list(loadings_samples = loadings_samples, cor_posterior = cor_posterior,
       corr_mean = corr_mean, corr_mean_zeroed = corr_mean_zeroed,
       prop_zeroed = mean(corr_mean_zeroed == 0))
}

results <- Map(function(f, YY) compute_basil_corr_zeroed(f, Y = YY, idx = idx),
               fits, Ys)

for (nm in names(results)) {
  dimnames(results[[nm]]$corr_mean)        <- list(colnames(Y_subset), colnames(Y_subset))
  dimnames(results[[nm]]$corr_mean_zeroed) <- list(colnames(Y_subset), colnames(Y_subset))
}
print(sapply(results, function(x) x$prop_zeroed))

# correlation heatmaps 
cols   <- scico(51, palette = "vik")
breaks <- seq(-1, 1, length.out = 52)

Corr_empirical <- cor(Y_subset)
emp_corr_plot  <- plot_correlation_heatmap(Corr_empirical, title = "Empirical",
                                           cluster = TRUE, breaks = breaks)
row_order <- emp_corr_plot$tree_row$order
col_order <- emp_corr_plot$tree_col$order

hm <- function(M, ttl) plot_correlation_heatmap(M, title = ttl, row_order = row_order,
                                                col_order = col_order, breaks = breaks)

LAB_INT <- "Normalized"

p_uq_orig <- hm(results[["Original"]]$corr_mean_zeroed, "Original")
p_uq_int  <- hm(results[["Rank-INT"]]$corr_mean_zeroed, LAB_INT)
p_pm_orig <- hm(results[["Original"]]$corr_mean,        "Original")
p_pm_int  <- hm(results[["Rank-INT"]]$corr_mean,        LAB_INT)

corrWB <- grid.arrange(p_uq_orig$gtable, p_uq_int$gtable,
                       p_pm_orig$gtable, p_pm_int$gtable, ncol = 2)
# ggsave('sensNorm_corr.jpeg', plot = corrWB, width = 9, height = 9)

## pathway gamma  
p_names <- colnames(allPaths)
idf_fa  <- 1:3

gamma_summary <- function(res) {
  G_s <- compute_Gamma_samples(res$loadings_samples$Lambda_samples, allPaths)
  lo  <- apply(G_s, c(1, 2), quantile, probs = 0.025)
  hi  <- apply(G_s, c(1, 2), quantile, probs = 0.975)
  G   <- apply(G_s, c(1, 2), mean)
  G[(lo < 0) & (hi > 0)] <- 0
  list(G = G, lo = lo, hi = hi, prop_zeroed = mean((lo < 0) & (hi > 0)))
}

gam <- lapply(results, gamma_summary)
print(sapply(gam, function(x) x$prop_zeroed))

TOP_N <- 5

build_union <- function(f) {
  top_idx <- function(G) order(abs(G[, f]), decreasing = TRUE)[seq_len(TOP_N)]
  idx_u   <- union(top_idx(gam$Original$G), top_idx(gam$`Rank-INT`$G))
  
  bind_rows(
    data.frame(Pathway = p_names[idx_u], Loading = gam$Original$G[idx_u, f],
               Factor_num = f, Data = "Original", stringsAsFactors = FALSE),
    data.frame(Pathway = p_names[idx_u], Loading = gam$`Rank-INT`$G[idx_u, f],
               Factor_num = f, Data = "Rank-INT", stringsAsFactors = FALSE)
  )
}

data_combined <- bind_rows(lapply(idf_fa, build_union)) %>%
  mutate(AbsLoading      = abs(Loading),
         Pathway_wrapped = str_wrap(str_trunc(gsub("_", " ", Pathway), 90), width = 30),
         Data            = factor(Data, levels = c("Original", "Rank-INT")))

# plot
make_factor_plot <- function(f_num) {
  df  <- data_combined %>% filter(Factor_num == f_num)
  ord <- df %>% group_by(Pathway_wrapped) %>%
    summarise(order_val = max(AbsLoading), .groups = "drop") %>%
    arrange(desc(order_val)) %>% pull(Pathway_wrapped)
  
  df <- df %>% mutate(rank  = match(Pathway_wrapped, ord),
                      y_num = rank + ifelse(Data == "Original", -0.15, 0.15))
  
  ggplot(df, aes(x = AbsLoading, y = y_num, shape = Data)) +
    geom_point(size = 2.5, colour = "black", stroke = 0.9) +
    scale_shape_manual(values = c("Original" = 16, "Rank-INT" = 4)) +
    scale_y_continuous(breaks = seq_along(ord), labels = ord, trans = "reverse") +
    coord_cartesian(xlim = c(0, max(data_combined$AbsLoading) * 1.05)) +
    labs(x = "Absolute loading", y = NULL, shape = NULL,
         title = paste0("Factor ", f_num)) +
    theme_minimal(base_size = 11) +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.position = "none",
          axis.text.y  = element_text(size = 8, lineheight = 0.8))
}

plots      <- lapply(idf_fa, make_factor_plot)
fplot_sens <- grid.arrange(grobs = plots, nrow = 1)

# ggsave('sens_dotplotNorm.jpeg', plot = fplot_sens, width = 11, height = 5)


