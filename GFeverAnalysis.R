library(BASIL)
library(tidyverse)
library(scico)
library(gridExtra)
library(PLIER)
library(dplyr)
library(tidyr)

rm(list=ls())

load("data_fever/data_Gfever.Rdata")
source("helper.R")
sc = colSums(geneSetMat)
pr1 = which(sc<11)
geneSetMat = geneSetMat[,-pr1]
sr = rowSums(geneSetMat)
gr1 = which(sr==0)
geneSetMat = geneSetMat[-gr1,]
data = data[-gr1,]

Y = t(data)
Y = scale(Y)

#-----------------# estimate BASIL #----------------------#
est_k = estimate_latent_dimension(Y, k_max=50)
fitGF = BASIL_point_estimates(Y,  geneSetMat, k=est_k$k_hat)
fitGF$tau_gamma/fitGF$tau_psi
Lambda_hat = fitGF$Lambda_C + fitGF$Lambda_N

#------------------------------------------------------------------------------#
#----------------------# correlation plot (Fig. 1) #---------------------------#

# select the 100 genes with largest variance
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


#------------------------------------------------------------------------------#
#----------------------------------# UQ for BASIL #----------------------------#

# UQ for Lambda
basil_loadings_samples = compute_posterior_samples_cc(
  Y, fitGF$Lambda_C, fitGF$Lambda_N, fitGF$tau_gamma, fitGF$tau_psi,
  fitGF$sigma_sq, fitGF$P_C, v0=1, sigma_sq_0=1, n_MC=200)

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
data_pplot = get_topPathways(Gamma_UQ[,1:10], p_names, top_n = 10, clean_names = TRUE)

# # rename labels too long
# lab_pth = data_pplot$Pathway
# lab_pth = sub(" ACTIVITY$", "", lab_pth)
# long_indices = which(nchar(lab_pth) > 30)
# lab_pth[long_indices]
# pathways_short = c(
#   "ABC GLUTATHIONE TRANSPORTER",              # 1
#   "PURINE NUCLEOTIDE TRANSPORTER",            # 2
#   "ATP INORGANIC ANION TRANS.",               # 3
#   "NUCLEOSIDE TRANSPORTER",                   # 4
#   "NUCLEOTIDE SUGAR TRANSPORTER",             # 5
#   "PYRIMIDINE NUCLEOSIDE TRANS.",             # 6
#   "ORGANOPHOSPHATE ANTIPORTER",               # 7
#   "L-GLUTAMATE TRANSPORTER",                  # 8
#   "PURINE NUCLEOTIDE TRANSPORTER",            # 9
#   "PYRIMIDINE NUCLEOSIDE TRANS.",             # 10
#   "L-GLUTAMATE TRANSPORTER",                  # 11
#   "ORGANOPHOSPHATE ANTIPORTER",               # 12
#   "ABC GLUTATHIONE TRANSPORTER",              # 13
#   "L-SERINE TRANSPORTER",                     # 14
#   "L-AMINO ACID TRANSPORTER",                 # 15
#   "L-LEUCINE TRANSPORTER",                    # 16
#   "NUCLEOSIDE TRANSPORTER",                   # 17
#   "VEGFR2 BINDING",                           # 18
#   "PURINE NUCLEOTIDE TRANSPORTER",            # 19
#   "VEGFR BINDING",                            # 20
#   "GABA-NA-CL SYMPORTER",                     # 21
#   "MONOATOMIC ANION TRANSPORTER",             # 22
#   "CHLORIDE TRANSPORTER",                     # 23
#   "ABC GLUTATHIONE TRANSPORTER",              # 24
#   "ATP INORGANIC ANION TRANS.",               # 25
#   "MONOATOMIC ANION TRANSPORTER",             # 26
#   "CHLORIDE TRANSPORTER",                     # 27
#   "HEME OXIDOREDUCTASE",                      # 28
#   "RETINOL DEHYDROG. (NADP+)",                # 29
#   "VEGFR BINDING",                            # 30
#   "VEGFR2 BINDING",                           # 31
#   "ATP CATION TRANSPORTER",                   # 32
#   "ACIDIC AMINO ACID TRANSPORTER",            # 33
#   "INTEGRIN BINDING (CELL-ECM)",              # 34
#   "PYRIMIDINE NUCLEOSIDE TRANS.",             # 35
#   "PURINE NUCLEOTIDE TRANSPORTER",            # 36
#   "NEUTRAL L-AA NA+ SYMPORTER",               # 37
#   "ATP INORGANIC ANION TRANS.",               # 38
#   "GLUTAMATE-GATED CA2+ CHANNEL",             # 39
#   "L-GLUTAMATE TRANSPORTER",                  # 40
#   "NUCLEOTIDE SUGAR TRANSPORTER",             # 41
#   "L-PROLINE TRANSPORTER",                    # 42
#   "PURINE NUCLEOTIDE TRANSPORTER",            # 43
#   "NUCLEOSIDE TRANSPORTER",                   # 44
#   "ABC GLUTATHIONE TRANSPORTER",              # 45
#   "ATP INORGANIC ANION TRANS.",               # 46
#   "NUCLEOTIDE SUGAR TRANSPORTER",             # 47
#   "ORGANOPHOSPHATE ANTIPORTER",               # 48
#   "ORGANIC HYDROXY TRANSPORTER",              # 49
#   "CARBOHYDRATE DERIV. TRANS.",               # 50
#   "L-GLUTAMATE TRANSPORTER",                  # 51
#   "POLYOL TRANSPORTER",                       # 52
#   "PURINE NUCLEOTIDE TRANSPORTER",            # 53
#   "ORGANOPHOSPHATE ESTER TRANS.",             # 54
#   "NUCLEOSIDE TRANSPORTER",                   # 55
#   "NUCLEOTIDE SUGAR TRANSPORTER",             # 56
#   "NUCLEOTIDE TRANSPORTER",                   # 57
#   "ABC GLUTATHIONE TRANSPORTER",              # 58
#   "ATP INORGANIC ANION TRANS.",               # 59
#   "CARBOHYDRATE DERIV. TRANS.",               # 60
#   "DEHYDROASCORBIC ACID TRANS.",              # 61
#   "PYRIMIDINE NUCLEOSIDE TRANS.",             # 62
#   "NUCLEOSIDE TRANSPORTER",                   # 63
#   "L-GLUTAMATE TRANSPORTER",                  # 64
#   "ABC GLUTATHIONE TRANSPORTER",              # 65
#   "PURINE NUCLEOTIDE TRANSPORTER",            # 66
#   "RETINOL DEHYDROG. (NADP+)",                # 67
#   "ATP INORGANIC ANION TRANS.",               # 68
#   "NUCLEOTIDE SUGAR TRANSPORTER",             # 69
#   "ORGANOPHOSPHATE ESTER TRANS.",             # 70
#   "ABC GLUTATHIONE TRANSPORTER",              # 71
#   "L-GLUTAMATE TRANSPORTER",                  # 72
#   "PYRIMIDINE NUCLEOSIDE TRANS.",             # 73
#   "NUCLEOSIDE TRANSPORTER",                   # 74
#   "ATP INORGANIC ANION TRANS.",               # 75
#   "ORGANOPHOSPHATE ANTIPORTER",               # 76
#   "NUCLEOTIDE SUGAR TRANSPORTER"              # 77
# )
# lab_pth[long_indices] = pathways_short
# max(nchar(lab_pth))
# data_pplot$Pathway = lab_pth

#--------# dotplot factor annontation #-----------#
fplot = dotplot_top_pathways(data_pplot)
fplot
# ggsave('Factor_annotGF.png', plot=fplot, device = 'png', width = 9.1, height = 9)

# posterior samples of Psi
Psi_samples = compute_Psi_samples(basil_loadings_samples$Lambda_samples, fitGF$P_C)

# 95% credible intervals
ci_lowerpsi = apply(Psi_samples, c(1, 2), quantile, probs = 0.025)
ci_upperpsi = apply(Psi_samples, c(1, 2), quantile, probs = 0.975)
contains_zeropsi = (ci_lowerpsi < 0) & (ci_upperpsi > 0)
table(contains_zeropsi)/length(contains_zeropsi) 
Psi_UQ = apply(Psi_samples, c(1, 2), mean)
Psi_UQ[contains_zero] = 0

normPsi = rep(0, nrow(Psi_UQ))
normPsi = apply(Psi_UQ, 1, function(x) sqrt(sum(x^2)))
top_Psi = order(normPsi, decreasing = TRUE)[1:20]
genePsi = rownames(geneSetMat[top_Psi,])
genePsi



#------------------------------------------------------------------------------#
#----------------------------------# Bar plot #--------------------------------#

c_sq = rowSums(fitGF$Lambda_C^2)[id]
p_sq = rowSums(fitGF$Lambda_N^2)[id]

gene_labs = colnames(Y)[id]


df = tibble(
  i = seq_along(c_sq),
  gene = gene_labs,
  CGamma = c_sq,
  Psi    = p_sq
) %>%
  mutate(total = CGamma + Psi,
         CGamma = CGamma / total,
         Psi    = Psi / total)  %>%
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
    legend.text = element_text(size = 25),
    plot.margin = margin(8, 10, 5, 20),
    axis.text.y = element_text(size = 15)
  )

var_explained_plot
#ggsave('var_explained_plot.png', plot=var_explained_plot, device = 'png', width = 15, height = 7.5)
