library(BASIL)
library(tidyverse)
library(scico)
library(gridExtra)
library(PLIER)
library(dplyr)
library(tidyr)

rm(list=ls())

# load("data_Gfever.Rdata")
source("functs/helper.R")

Y = t(data)
Y = scale(Y)

#-----------------# estimate BASIL #----------------------#
est_k = estimate_latent_dimension(Y, k_max=50)
fitGF = BASIL_point_estimates(Y,  geneSetMat, k=est_k$k_hat)
fitGF$tau_gamma/fitGF$tau_psi
Lambda_hat = fitGF$Lambda_C + fitGF$Lambda_N

#------------------------------------------------------------------------------#
#----------------------# correlation plot (Fig. 1) #---------------------------#

# subset 100 genes
set.seed(134)
id <- sample(c(1:dim(Y)[2]), 100)
sub_data <- Y[,id]

# empirical correlation
corP = cor(sub_data)

# BASIl correlation
CovBasil = Lambda_hat %*% t(Lambda_hat) + fitGF$sigma_sq*diag(ncol(Y))
CovBasil_sub = CovBasil[id,id]
corrBasil = cov2cor(CovBasil_sub)

# fit PLIER
PLresult = PLIER(t(Y), geneSetMat, k = est_k$k_hat, doCrossval = F, scale = F)
covPLIER = PLresult$Z %*% (cov(t(PLresult$B))) %*% t(PLresult$Z)
sigma_sqPLIER = mean((Y - t(PLresult$B) %*% t(PLresult$Z))^2)
covPLIER = covPLIER + sigma_sqPLIER * diag(ncol(Y))
CovPLIER_sub = covPLIER[id,id]
corrPLIER = cov2cor(CovPLIER_sub)

#--------# correlation plot #-------#
cols = scico(51, palette = "vik") # color blind friendly palette
breaks = seq(-1, 1, length.out = 52)

# empirical correlation
emp_corr_plot = plot_correlation_heatmap(corP, title = "Empirical", cluster = TRUE,
                                         breaks = breaks)
row_order = emp_corr_plot$tree_row$order
col_order = emp_corr_plot$tree_col$order

# BASIL
basil_corr_plot = plot_correlation_heatmap(corrBasil, title = "BASIL", row_order = row_order, 
  col_order = col_order, breaks = breaks)
# PLIER
plier_corr_plot = plot_correlation_heatmap(corrPLIER, title = "PLIER", row_order = row_order, 
  col_order = col_order, breaks = breaks)
corrNP = grid.arrange(emp_corr_plot$gtable, basil_corr_plot$gtable,
                               plier_corr_plot$gtable, ncol = 3)

# ggsave('CorrFever.png', plot=corrNP, device = 'png', width = 9.1, height = 3.3)

#------# observed vs predicted plot #--------#

df_basil_cplot = data.frame(Observed = corP[lower.tri(corP)],
                            Predicted = corrBasil[lower.tri(corrBasil)],
                            Method = "BASIL")
df_plier_cplot = data.frame(Observed = corP[lower.tri(corP)],
                            Predicted = corrPLIER[lower.tri(corrPLIER)],
                            Method = "PLIER")
lims = range(c(df_basil_cplot$Observed, df_basil_cplot$Predicted,
               df_plier_cplot$Observed, df_plier_cplot$Predicted))

p1_obs = plot_correlation_scatter(df_basil_cplot, title = "BASIL", lims)
p2_obs = plot_correlation_scatter(df_plier_cplot, title = "PLIER", lims)

obspPLot = grid.arrange(p1_obs, p2_obs, nrow=1)
# ggsave('ObsPFever.png', plot=obspPLot, device = 'png', width = 8, height = 4)



#------------------------------------------------------------------------------#
#----------------------------------# UQ for BASIL #----------------------------#

# UQ for Lambda
basil_loadings_samples = compute_posterior_samples_cc(
  Y, fitGF$Lambda_C, fitGF$Lambda_N, fitGF$tau_gamma, fitGF$tau_psi,
  fitGF$sigma_sq, fitGF$P_C, v0=1, sigma_sq_0=1, n_MC=500)
# save(basil_loadings_samples, file="UQ_BASIL.Rdata")

# posterior samples of Gamma
Gamma_samples = compute_Gamma_samples(basil_loadings_samples$Lambda_samples,
                                          geneSetMat)

# 95% credible intervals
ci_lower = apply(Gamma_samples, c(1, 2), quantile, probs = 0.025)
ci_upper = apply(Gamma_samples, c(1, 2), quantile, probs = 0.975)
contains_zero = (ci_lower < 0) & (ci_upper > 0)
table(contains_zero)/length(contains_zero)
Gamma_UQ = apply(Gamma_samples, c(1, 2), mean)
Gamma_UQ[contains_zero] = 0

# get top pathways
p_names = colnames(geneSetMat)
data_pplot = get_topPathways(Gamma_UQ[,1:10], p_names, top_n = 5, 
                             clean_names = TRUE, ci_lower = ci_lower[,1:10],
                             ci_upper = ci_upper[,1:10])

#--------# dotplot factor annontation #-----------#
lab_pth = data_pplot$Pathway
lab_pth = sub(" ACTIVITY$", "", lab_pth)
long_indices <- which(nchar(lab_pth) > 30)
lab_pth[long_indices]
shortened_names  <- c(
  "PI(3,5)P2 3-PHOSPHATASE",         # 1
  "PI BISPHOSPHATE PHOSPHATASE",     # 2
  "PI PHOSPHATE PHOSPHATASE",        # 3
  "PI MONOPHOSPHATE PHOSPHATASE",    # 4
  "PI(3,5)P2 5-PHOSPHATASE",         # 5
  "MACROMOLECULE TRANSPORTER",       # 6
  "ABC GLUTATHIONE TRANSPORTER",     # 7
  "PI MONOPHOSPHATE PHOSPHATASE",    # 8
  "PI(3,5)P2 3-PHOSPHATASE",         # 9
  "PROTEIN TRANSPORTER",             # 10
  "POSTSYNAPTIC NT RECEPTOR",        # 11
  "DNA DAMAGE SENSOR",               # 12
  "BETA-2 ADRENERGIC BINDING",       # 13
  "ANION ATPASE",                    # 14
  "ABC GLUTATHIONE TRANSPORTER",     # 15
  "CATION ATPASE",                   # 16
  "CALCIUM ANTIPORTER",              # 17
  "PI PHOSPHATE PHOSPHATASE",        # 18
  "PI(3,5)P2 3-PHOSPHATASE",         # 19
  "NAD(P)H HEME OXIDOREDUCTASE",     # 20
  "PYRIMIDINE TRANSPORT",            # 21
  "PI(3,5)P2 3-PHOSPHATASE",         # 22
  "PI MONOPHOSPHATE PHOSPHATASE",    # 23
  "PI PHOSPHATE PHOSPHATASE",        # 24
  "PI(3,5)P2 PHOSPHATASE",           # 25
  "PKA REGULATOR",                   # 26
  "PI(3,5)P2 3-PHOSPHATASE",         # 27
  "PI MONOPHOSPHATE PHOSPHATASE",    # 28
  "PI PHOSPHATE PHOSPHATASE",        # 29
  "PI(3,5)P2 5-PHOSPHATASE",         # 30
  "CATION ATPASE",                   # 31
  "ABC GLUTATHIONE TRANSPORTER",     # 32
  "PI(3,5)P2 3-PHOSPHATASE",         # 33
  "PI(3,5)P2 PHOSPHATASE",           # 34
  "PI MONOPHOSPHATE PHOSPHATASE",    # 35
  "MACROMOLECULE TRANSPORTER",       # 36
  "PI(3,5)P2 5-PHOSPHATASE",         # 37
  "PI(3,5)P2 3-PHOSPHATASE",         # 38
  "POSTSYNAPTIC NT RECEPTOR",        # 39
  "PI MONOPHOSPHATE PHOSPHATASE",    # 40
  "PI PHOSPHATE PHOSPHATASE",        # 41
  "PI(4,5)P2 5-PHOSPHATASE"          # 42
)
  
lab_pth[long_indices] = shortened_names
max(nchar(lab_pth))
data_pplot$Pathway = lab_pth

fplot = dotplot_top_pathways(data_pplot)
fplot
# ggsave('Factor_annotGF_new.png', plot=fplot, device = 'png', width = 9.1, height = 9)


#------------------------------------------------------------------------------#
#----------------------------------# Bar plot #--------------------------------#

gene_variances <- apply(t(data), 2, var)
id50 <- order(gene_variances, decreasing = TRUE)[1:50]
sub_data <- Y[,id50]

c_sq = rowSums(fitGF$Lambda_C^2)[id50]
p_sq = rowSums(fitGF$Lambda_N^2)[id50]
gene_labs = colnames(Y)[id50]


df = tibble(i = seq_along(c_sq), gene = gene_labs, CGamma = c_sq, Psi = p_sq) %>%
  mutate(total = CGamma + Psi, CGamma = CGamma / total, Psi    = Psi / total)  %>%
  pivot_longer(cols = c(CGamma, Psi), names_to = "component", values_to = "prop")

df = df %>%
  group_by(gene) %>%
  arrange(desc(component == "CGamma"), .by_group = TRUE) %>%
  ungroup()

df$gene = factor(
  df$gene,
  levels = df %>%
    filter(component == "CGamma") %>%
    arrange(desc(prop)) %>%
    pull(gene)
)
df$component <- factor(df$component, levels = c("CGamma", "Psi"))

var_explained_plot <- ggplot(df, aes(x = gene, y = prop, fill = component)) +
  geom_col(width = 0.9, position = position_stack(reverse = TRUE)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0001)) +
  scale_fill_manual(
    name = "",
    values = c(CGamma = "#6A3D9A", Psi = "#33A02C"),
    labels = c(expression(C~Gamma), expression(Psi))
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
    legend.text = element_text(size = 30),
    plot.margin = margin(8, 10, 5, 20),
    axis.text.y = element_text(size = 15)
  )

var_explained_plot

# ggsave('var_explained_plot50.png', plot=var_explained_plot, device = 'png', 
#        width = 15, height = 7.5)


#--------------------------------# network plot #------------------------------#

# idx <- order(gene_variances, decreasing = TRUE)[1:100]
# CovPLIER_subx = covPLIER[idx,idx]
# corrPLIERx = cov2cor(CovPLIER_subx)
# 
# # Compute correlation posterior 
# cor_posterior <- compute_correlation_posterior_samples_cc(
#   Lambda_samples = basil_loadings_samples$Lambda_samples[idx,,],
#   sigma_sq_samples = basil_loadings_samples$sigma_sq_samples,
#   samples = TRUE
# )
# CorrBasil_mean <- cor_posterior$posterior_mean
# CorrBasil_mean_zeroed <- CorrBasil_mean
# # Compute 95% credible intervals
# corrBasil_qs <- apply(cor_posterior$posterior_samples, c(1, 2), quantile, 
#                       probs = c(0.025, 0.975))
# CorrBasil_mean_zeroed[((corrBasil_qs[1,,] <0) & (corrBasil_qs[2,,] >0))] = 0
# # Proportion of correlations set to zero
# mean(CorrBasil_mean_zeroed==0)
# 
# Y_subset <- Y[,idx]
# row.names(CorrBasil_mean_zeroed) <- colnames(Y_subset)
# colnames(CorrBasil_mean_zeroed) <- colnames(Y_subset)
# 
# network_plotBASIL <- plot_gene_network(CorrBasil_mean_zeroed, n = 50, max_overlaps = 50)
# network_plotBASIL
# 
# 
# network_plotPLIER <- plot_gene_network(corrPLIERx, n = 50)
# network_plotPLIER
