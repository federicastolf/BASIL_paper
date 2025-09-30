library(tidyverse)
library(scico)
library(gridExtra)

rm(list=ls())

load("data_fever/data_Gfever.Rdata")
source('helper_functions.R')

Y = t(data)
Y = scale(Y)

est_k = estimate_latent_dimension(Y, k_max=50)
fitGF = compute_point_estimates(Y,  geneSetMat, k=est_k$k_hat)
fitGF$tau_C/fitGF$tau_N

Lambda_hat = fitGF$Lambda_C + fitGF$Lambda_N
CovBasil = Lambda_hat %*% t(Lambda_hat) + fitGF$sigma_sq*diag(ncol(Y))
corrBasil = cov2cor(CovBasil)

#----------------# check grafico correlation #--------------------#
# random select 50 genes
set.seed(139)
id = sample(c(1:NCOL(Y)), 50, replace = F)
sub_data = Y[,id]

# empirical correlation 
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

sub_corrBasil = corrBasil[id,id]
BSPlot = pheatmap(sub_corrBasil[row_order, col_order], treeheight_row = 0, treeheight_col = 0,
                  cluster_rows = F, cluster_cols = F, border_color ="NA", legend=F,
                  show_colnames = F, show_rownames = F, breaks = breaks, color = cols, 
                  fontsize = 12, main = "BASIL")

CorrNP = grid.arrange(EmpPlot$gtable, BSPlot$gtable,  ncol = 2)
