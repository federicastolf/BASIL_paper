library(PLIER)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(pheatmap)
library(scico)
library(mvtnorm)
rm(list=ls())

#---------------------------# human whole blood data #------------------------#

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
dataWholeBlood_top100 <- t(dataWholeBlood)[,id]


#----# fit PLIER #-----#
PLresult = PLIER(t(Y), allPaths, k = 8, doCrossval = F, scale = F)
covPLIER = PLresult$Z %*% (cov(t(PLresult$B))) %*% t(PLresult$Z)
corrPLIER = cov2cor(covPLIER)

#----# fit BASIL #-----#
source('helper_functions.R')
svd_Y = svd(Y)
k_max = min(which((cumsum(svd_Y$d^2) / sum(svd_Y$d^2) > 0.8)))
k_max
est_kWB = estimate_latent_dimension(Y, k_max=k_max)
fitWB = compute_point_estimates(Y,  allPaths, k=est_kWB$k_hat)
fitWB$tau_C/fitWB$tau_N

Lambda_hatWB = fitWB$Lambda_C + fitWB$Lambda_N
CovBasilWB = Lambda_hatWB %*% t(Lambda_hatWB) + fitWB$sigma_sq*diag(ncol(Y))
corrBasilWb = cov2cor(CovBasilWB)

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

#------------------------# correlation Plot #---------------------#

# empirical correlation 
sub_data <- Y[,id]
corP = cor(sub_data)

# color blind palette
cols = scico(51, palette = "vik") 
breaks = seq(-1, 1, length.out = 52)
EmpPlot = pheatmap(corP, treeheight_row = 0, treeheight_col = 0, cluster_rows = T,
                   cluster_cols = T, border_color ="NA", legend=F, show_colnames = F,
                   show_rownames = F, breaks = breaks, color = cols, fontsize = 12,
                   main = "Empirical")

row_order = EmpPlot$tree_row$order
col_order = EmpPlot$tree_col$order

row.names(corrBasil_mean) = colnames(sub_data)
colnames(corrBasil_mean) = colnames(sub_data)


BSPlot = pheatmap(corrBasil_mean[row_order, col_order], treeheight_row = 0, treeheight_col = 0,
                  cluster_rows = F, cluster_cols = F, border_color ="NA", legend=F,
                  show_colnames = F, show_rownames = F, breaks = breaks, color = cols, 
                  fontsize = 12, main = "BASIL")

row.names(corrBasil_mean_zeroed) = colnames(sub_data)
colnames(corrBasil_mean_zeroed) = colnames(sub_data)

BSPlot_uq = pheatmap(corrBasil_mean_zeroed[row_order, col_order], treeheight_row = 0, treeheight_col = 0,
                     cluster_rows = F, cluster_cols = F, border_color ="NA", legend=F,
                     show_colnames = F, show_rownames = F, breaks = breaks, color = cols, 
                     fontsize = 12, main = "BASIL-UQ")


sub_corrPLIER = corrPLIER[id,id]
PLPlot = pheatmap(sub_corrPLIER[row_order, col_order], treeheight_row = 0, treeheight_col = 0, 
                  cluster_rows = F, cluster_cols = F, border_color ="NA", legend=F, 
                  show_colnames = F, show_rownames = F, breaks = breaks, color = cols, 
                  fontsize = 12, main="PLIER")

grid.arrange(EmpPlot$gtable, BSPlot$gtable, BSPlot_uq$gtable,  PLPlot$gtable, ncol = 4)
CorrNP = grid.arrange(EmpPlot$gtable, BSPlot$gtable, BSPlot_uq$gtable,  PLPlot$gtable, ncol = 4)
# ggsave('CorrWholeBLOOD1.png', plot=CorrNP, device = 'png', width = 9.1, height = 3.3)


#-------------------# observed vs predicted plot #-------------------------#


df_basil = data.frame(Observed = corP[lower.tri(corP)], 
                      Predicted = corrBasil_mean[lower.tri(corrBasil_mean)], 
                      Method = "BASIL")
df_plier = data.frame(Observed = corP[lower.tri(corP)], 
                      Predicted = sub_corrPLIER[lower.tri(sub_corrPLIER)],
                      Method = "PLIER")
lims = range(c(df_basil$Observed, df_basil$Predicted, 
               df_plier$Observed, df_plier$Predicted))

# library(scales)
# library(ggthemes)
# show_col(tableau_color_pal("Color Blind")(10))

# BASIL plot
p1 = ggplot(df_basil, aes(x = Observed, y = Predicted)) +
  geom_point(color = "#1170aa", shape = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1) +
  labs(title = "BASIL", x = "Observed", y = "Predicted") +
  xlim(lims) + ylim(lims) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.title=element_text(size=16))
p1

# PLIER plot
p2 = ggplot(df_plier, aes(x = Observed, y = Predicted)) +
  geom_point(color = "#1170aa", shape = 1) +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1) +
  labs(title = "PLIER", x = "Observed", y = "Predicted") +
  xlim(lims) + ylim(lims) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.title=element_text(size=16),)
p2

obspPLot = grid.arrange(p1,p2, nrow=1)

dev.new()
obspPLot

# ggsave('ObsPWholeBLOOD.png', plot=obspPLot, device = 'png', width = 8, height = 4)


#-------------------# gene network plot #-------------------------#

library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(reshape2)


df_long <- melt(corrBasil_mean[]) %>%
  filter(value < 0.9999)
df_sorted <- df_long %>%
  arrange(desc(value))
head(df_sorted, n=20)


dev.new()
network_plot <- plot_gene_network(corrBasil_mean_zeroed, n = 50)
ggsave('network_plot_WholeBLOOD.png', plot=network_plot, device = 'png', width = 8, height = 4)




