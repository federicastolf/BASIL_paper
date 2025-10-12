library(tidyverse)
library(scico)
library(gridExtra)
library(pheatmap)
library(PLIER)

rm(list=ls())

load("data_fever/data_Gfever.Rdata")
source('helper_functions.R')
#load("data_fever/fever_covariate.Rdata")
source("plot_functions.R")

Y = t(data)
Y = scale(Y)

#-----------------# estimate BASIL #----------------------#
est_k = estimate_latent_dimension(Y, k_max=50)
fitGF = compute_point_estimates(Y,  geneSetMat, k=est_k$k_hat)
fitGF$tau_C/fitGF$tau_N
Lambda_hat = fitGF$Lambda_C + fitGF$Lambda_N

#------------------------------------------------------------------------------#
#----------------------# correlation plot (Fig. 1) #---------------------------#

# selet the 100 genes with largest variance
gene_variances <- apply(Y, 2, var)
id <- order(gene_variances, decreasing = TRUE)[1:100]
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
CovPLIER_sub = covPLIER[id,id]
corrPLIER = cov2cor(CovPLIER_sub)

#--------# correlation plot #-------#
cols = scico(51, palette = "vik") # color blind friendly palette
breaks = seq(-1, 1, length.out = 52)

emp_corr_plot = plot_correlation_heatmap(corP, title = "Empirical", cluster = TRUE)
row_order = emp_corr_plot$tree_row$order
col_order = emp_corr_plot$tree_col$order
basil_corr_plot = plot_correlation_heatmap(corrBasil, title = "BASIL", 
                                           row_order = row_order, col_order = col_order)
plier_corr_plot = plot_correlation_heatmap(corrPLIER, title = "PLIER",
                                           row_order = row_order, col_order = col_order)
corrNP = grid.arrange(emp_corr_plot$gtable, basil_corr_plot$gtable, 
                               plier_corr_plot$gtable, ncol = 3)
# ggsave('CorrFever.png', plot=CorrNP, device = 'png', width = 9.1, height = 3.3)

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

