
library(PLIER)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(pheatmap)
library(scico)
rm(list=ls())

source('helper_functions.R')
source("plot_functions.R")

# load data
data("dataWholeBlood") 
data("bloodCellMarkersIRISDMAP")
data("canonicalPathways")
allPaths = combinePaths(bloodCellMarkersIRISDMAP, canonicalPathways)
cm = intersect(rownames(dataWholeBlood), rownames(allPaths))
allPaths = allPaths[cm,]
dataWholeBlood=dataWholeBlood[cm,]
Y = t(dataWholeBlood)
Y = scale(Y)


#------------------------------------------------------------------------------#
#--------------------------# Correlation plot (Fig.1) #------------------------#

# fit PLIER 
PL10 = PLIER(t(Y), allPaths, k = 10, doCrossval = F, scale = F)
covPLIER = PL10$Z %*% (cov(t(PL10$B))) %*% t(PL10$Z)
corrPLIER = cov2cor(covPLIER)
# fit BASIL 
BL10 = compute_point_estimates(Y,  allPaths, k=10)
Lambda_hat = BL10$Lambda_C + BL10$Lambda_N
CovBasil = Lambda_hat %*% t(Lambda_hat) + BL10$sigma_sq*diag(ncol(Y))
corrBasil = cov2cor(CovBasil)

# random select 50 genes
set.seed(139)
id = sample(c(1:NCOL(Y)), 50, replace = F)
sub_data = Y[,id]
# empirical correlation 
corP = cor(sub_data)
sub_corrBasil = corrBasil[id,id]
sub_corrPLIER = corrPLIER[id,id]

#--------# correlation plot #-------#
cols = scico(51, palette = "vik") # color blind friendly palette
breaks = seq(-1, 1, length.out = 52)

emp_corr_plot = plot_correlation_heatmap(corP, title = "Empirical", cluster = TRUE)
row_order = emp_corr_plot$tree_row$order
col_order = emp_corr_plot$tree_col$order
basil_corr_plot = plot_correlation_heatmap(sub_corrBasil, title = "BASIL", 
                                      row_order = row_order, col_order = col_order)
plier_corr_plot = plot_correlation_heatmap(sub_corrPLIER, title = "PLIER",
                                       row_order = row_order, col_order = col_order)

corr_comparison = grid.arrange(emp_corr_plot$gtable, basil_corr_plot$gtable, 
                                plier_corr_plot$gtable, ncol = 3)
# ggsave('CorrWholeBLOOD1.png', plot=CorrNP, device = 'png', width = 9.1, height = 3.3)

#------# observed vs predicted plot #--------#

df_basil_cplot = data.frame(Observed = corP[lower.tri(corP)], 
                      Predicted = sub_corrBasil[lower.tri(sub_corrBasil)], 
                      Method = "BASIL")
df_plier_cplot = data.frame(Observed = corP[lower.tri(corP)], 
                       Predicted = sub_corrPLIER[lower.tri(sub_corrPLIER)],
                       Method = "PLIER")
lims = range(c(df_basil_cplot$Observed, df_basil_cplot$Predicted,
                df_plier_cplot$Observed, df_plier_cplot$Predicted))

p1_obs = plot_correlation_scatter(df_basil_cplot, title = "BASIL", lims)
p2_obs = plot_correlation_scatter(df_plier_cplot, title = "PLIER", lims)

obspPLot = grid.arrange(p1_obs, p2_obs, nrow=1)

# ggsave('ObsPWholeBLOOD.png', plot=obspPLot, device = 'png', width = 8, height = 4)

