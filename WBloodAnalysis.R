library(PLIER)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(pheatmap)
library(scico)
library(mvtnorm)
library(circlize)
library(ComplexHeatmap)
library(ggraph)
library(igraph)
library(RColorBrewer)

rm(list=ls())

source('helper_functions.R')
source('plot_functions.R')

# load data
data("dataWholeBlood") 
data("bloodCellMarkersIRISDMAP")
data("canonicalPathways")
allPaths = combinePaths(bloodCellMarkersIRISDMAP, canonicalPathways)
cm=intersect(rownames(dataWholeBlood), rownames(allPaths))
allPaths=allPaths[cm,]
dataWholeBlood=dataWholeBlood[cm,]
Y = t(dataWholeBlood)
Y = scale(Y, scale=T)
n = nrow(Y)
p = ncol(Y)

# filter 100 genes with largest variance
gene_variances <- apply(t(dataWholeBlood), 2, var)
id <- order(gene_variances, decreasing = TRUE)[1:100]

#----# fit PLIER #-----#
PLresult = PLIER(t(Y), allPaths, k = 8, doCrossval = F, scale = F)
covPLIER = PLresult$Z %*% (cov(t(PLresult$B))) %*% t(PLresult$Z)
corrPLIER = cov2cor(covPLIER)

#----# fit BASIL #-----#
svd_Y = svd(Y)
k_max = min(which((cumsum(svd_Y$d^2) / sum(svd_Y$d^2) > 0.8)))
k_max
est_kWB = estimate_latent_dimension(Y, k_max=k_max)
fitWB = compute_point_estimates(Y,  allPaths, k=est_kWB$k_hat)
fitWB$tau_C/fitWB$tau_N

Lambda_hatWB = fitWB$Lambda_C + fitWB$Lambda_N
CovBasilWB = Lambda_hatWB %*% t(Lambda_hatWB) + fitWB$sigma_sq*diag(ncol(Y))
corrBasilWb = cov2cor(CovBasilWB)

# UQ for BASIL
basil_loadings_samples = compute_posterior_samples_cc(
  Y, fitWB$Lambda_C, fitWB$Lambda_N, fitWB$tau_C, fitWB$tau_N,
  fitWB$sigma_sq, fitWB$P_C, v0=1, sigma_sq_0=1, n_MC=500
)

cor_est = compute_correlation_posterior_samples_cc(
  basil_loadings_samples$Lambda_samples[id,,], basil_loadings_samples$sigma_sq_samples, samples=T
)

corrBasil_mean = cor_est$posterior_mean
corrBasil_mean_zeroed = corrBasil_mean
corrBasil_qs = apply(cor_est$posterior_samples, c(1,2), function(x)(quantile(x, probs=c(0.025, 0.975))))
corrBasil_mean_zeroed[((corrBasil_qs[1,,] <0) & (corrBasil_qs[2,,] >0))] = 0
mean(corrBasil_mean_zeroed==0)

#---------------------------------# correlation Plot #-------------------------#

# empirical correlation 
sub_data <- Y[,id]
corP = cor(sub_data)

# color blind palette
cols = scico(51, palette = "vik") 
breaks = seq(-1, 1, length.out = 52)

# empirical correlation
emp_corr_plot = plot_correlation_heatmap(corP, title = "Empirical", cluster = TRUE)
row_order = emp_corr_plot$tree_row$order
col_order = emp_corr_plot$tree_col$order

# BASIL
row.names(corrBasil_mean) = colnames(sub_data)
colnames(corrBasil_mean) = colnames(sub_data)
basil_corr_plot = plot_correlation_heatmap(corrBasil_mean, title = "BASIL", 
                                           row_order = row_order, col_order = col_order)

# BASIL UQ
row.names(corrBasil_mean_zeroed) = colnames(sub_data)
colnames(corrBasil_mean_zeroed) = colnames(sub_data)
basilUQ_corr_plot = plot_correlation_heatmap(
  corrBasil_mean_zeroed, title = "BASIL-UQ", row_order = row_order, col_order = col_order)

#PLIER
sub_corrPLIER = corrPLIER[id,id]
plier_corr_plot = plot_correlation_heatmap(sub_corrPLIER, title = "PLIER",
                                           row_order = row_order, col_order = col_order)

corrWB = grid.arrange(emp_corr_plot$gtable, basil_corr_plot$gtable, basilUQ_corr_plot$gtable,
                               plier_corr_plot$gtable, ncol = 4)
# ggsave('CorrWholeBLOOD1.png', plot=CorrNP, device = 'png', width = 9.1, height = 3.3)

#----------------------------# gene network plot #-----------------------------#

df_long <- melt(corrBasil_mean[]) %>%
  filter(value < 0.9999)
df_sorted <- df_long %>%
  arrange(desc(value))
head(df_sorted, n=20)


network_plot <- plot_gene_network(corrBasil_mean_zeroed, n = 50)
network_plot
# ggsave('network_plot_WholeBLOOD.png', plot=network_plot, device = 'png', 
#        width = 8, height = 4)




