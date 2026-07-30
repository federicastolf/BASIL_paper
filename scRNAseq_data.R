
# BiocManager::install("scRNAseq")
library(scRNAseq)
library(scater)
library(msigdbr)
library(BASIL)
library(scico)
library(dplyr)
library(tidyr)
library(gridExtra)

rm(list=ls())

#----------------------------# download data #--------------------------------#

sce = ZeiselBrainData()

# merge duplicate rows (same gene, different genomic loci -> "_loc" suffix)
sce = aggregateAcrossFeatures(sce, id = sub("_loc[0-9]+$", "", rownames(sce)))
counts_mat = counts(sce)   # genes x cells
dim(counts_mat)

# get REACTOME pathways for mouse 
reactome_sets = msigdbr(species = "Mus musculus", category = "C2", 
                         subcategory = "CP:REACTOME")

# keep only gene-pathway pairs for genes actually present in the data
reactome_sets = reactome_sets[reactome_sets$gene_symbol %in% rownames(counts_mat), ]
C = table(reactome_sets$gene_symbol, reactome_sets$gs_name)
C = as.matrix(C)
dim(C)   # genes x pathways

common_genes = intersect(rownames(counts_mat), rownames(C))
counts_mat = counts_mat[common_genes, ]
C = C[common_genes, ]
dim(counts_mat)   # genes x cells
dim(C)      # genes x pathways

# filtered out pathways with fewer than 10 genes and genes that were not
# annotated to any gene set
sc = colSums(C)
pr1 = which(sc<11)
C = C[,-pr1]
sr = rowSums(C)
gr1 = which(sr==0)
C = C[-gr1,]
counts_mat = counts_mat[-gr1,]


#-----# data standardization #-----#

Ys = t(counts_mat)
Ys = scale(log1p(Ys))
dim(Ys)
summary(c(Ys))

#----------------------------# fit BASIL #-------------------------------------#

est_k = estimate_latent_dimension(Ys, k_max=50) # 19
fitsce = BASIL_point_estimates(Ys,  C, k=est_k$k_hat)
fitsce$tau_gamma/fitsce$tau_psi
Lambda_hat = fitsce$Lambda_C + fitsce$Lambda_N

# UQ for Lambda
basil_loadings_samples = compute_posterior_samples_cc(
  Ys, fitsce$Lambda_C, fitsce$Lambda_N, fitsce$tau_gamma, fitsce$tau_psi,
  fitsce$sigma_sq, fitsce$P_C, v0=1, sigma_sq_0=1, n_MC=500)

#---------# correlation plot #------------#
gene_variances = apply(t(counts_mat), 2, var)
idx = order(gene_variances, decreasing = TRUE)[1:100]

cols = scico(51, palette = "vik") 
breaks = seq(-1, 1, length.out = 52)

# empirical correlation
sub_data = t(counts_mat)[,idx]
corP = cor(sub_data)
emp_corr_plot = plot_correlation_heatmap(corP, title = "Empirical", cluster = TRUE,
                                         breaks = breaks)
row_order = emp_corr_plot$tree_row$order
col_order = emp_corr_plot$tree_col$order

# BASIL correlation
cor_posterior = compute_correlation_posterior_samples_cc(
  Lambda_samples = basil_loadings_samples$Lambda_samples[idx,,],
  sigma_sq_samples = basil_loadings_samples$sigma_sq_samples,
  samples = TRUE
)

CorrBasil_mean = cor_posterior$posterior_mean
row.names(CorrBasil_mean) = colnames(sub_data)
colnames(CorrBasil_mean) = colnames(sub_data)
basil_corr_plot = plot_correlation_heatmap(
  CorrBasil_mean, title = "BASIL", row_order = row_order, col_order = col_order, 
  breaks = breaks
)

corr_sce = grid.arrange(emp_corr_plot$gtable, basil_corr_plot$gtable, ncol = 2)
# ggsave('CorrSce.png', plot=corr_sce, device = 'png', width = 9, height = 4.7)


#----------# factor annotation #---------#

# posterior samples of Gamma
Gamma_samples = compute_Gamma_samples(basil_loadings_samples$Lambda_samples, C)

# 95% credible intervals
ci_lower = apply(Gamma_samples, c(1, 2), quantile, probs = 0.025)
ci_upper = apply(Gamma_samples, c(1, 2), quantile, probs = 0.975)
contains_zero = (ci_lower < 0) & (ci_upper > 0)
table(contains_zero)/length(contains_zero)
Gamma_UQ = apply(Gamma_samples, c(1, 2), mean)
Gamma_UQ[contains_zero] = 0

# get top pathways
p_names = colnames(C)
idf_fa = c(1:10)
data_pplot = get_topPathways(Gamma_UQ[,idf_fa], p_names, top_n = 5, 
                             clean_names = TRUE, ci_lower = ci_lower[,idf_fa],
                             ci_upper = ci_upper[,idf_fa])

fplot = dotplot_top_pathways(data_pplot)

lab_pth = data_pplot$Pathway
long_indices = which(nchar(lab_pth) > 30)
lab_pth[long_indices]

shortened_names = c(
  "TBK1/IKK-IRF3/IRF7 ACTIVATION",    # 1
  "TICAM1-IRF3/IRF7 ACTIVATION",      # 2
  "BCR SIGNALING",                    # 3
  "BCR ANTIGEN 2ND MESSENGERS",       # 4
  "GOLGI-ER RETROGRADE TRANSPORT",    # 5
  "EGFR CANCER VARIANT SIGNALING",    # 6
  "EGFRvIII CONSTITUTIVE SIGNALING",  # 7
  "EGFR CANCER VARIANT SIGNALING",    # 8
  "BCR DOWNSTREAM SIGNALING",         # 9
  "EGFRvIII CONSTITUTIVE SIGNALING",  # 10
  "BCR ANTIGEN 2ND MESSENGERS",       # 11
  "BCR SIGNALING",                    # 12
  "TICAM1-IRF3/IRF7 ACTIVATION",      # 13
  "TBK1/IKK-IRF3/IRF7 ACTIVATION",    # 14
  "RIP1-MEDIATED IKK RECRUITMENT",    # 15
  "FGFR4 DOWNSTREAM SIGNALING",       # 16
  "EGFR CANCER VARIANT SIGNALING",    # 17
  "FGFR2 DOWNSTREAM SIGNALING",       # 18
  "SPRY REGULATION OF FGF",           # 19
  "EGFR CANCER VARIANT SIGNALING",    # 20
  "EGFRvIII CONSTITUTIVE SIGNALING",  # 21
  "PD-L1 EXPRESSION REG",             # 22
  "PD-L1 POST-TRANSLATIONAL REG",     # 23
  "SPRY REGULATION OF FGF",           # 24
  "EGFRvIII CONSTITUTIVE SIGNALING",  # 25
  "EGFR CANCER VARIANT SIGNALING",    # 26
  "ADP SIGNALING VIA P2Y12",          # 27
  "TICAM1-RIP1 IKK RECRUITMENT",      # 28
  "RIP1-MEDIATED IKK RECRUITMENT",    # 29
  "TBK1/IKK-IRF3/IRF7 ACTIVATION",    # 30
  "TICAM1-IRF3/IRF7 ACTIVATION",      # 31
  "TELOMERE LAGGING STRAND SYNTHESIS",# 32
  "NOTCH3 ICD TRANSCRIPTION",         # 33
  "NOTCH3 SIGNAL TO NUCLEUS",         # 34
  "CDC25A UBIQ. DEGRADATION",         # 35
  "TBK1/IKK-IRF3/IRF7 ACTIVATION",    # 36
  "RIP1-MEDIATED IKK RECRUITMENT",    # 37
  "P53-DEPENDENT G1 DNA DAMAGE"       # 38
)

lab_pth[long_indices] = shortened_names
max(nchar(lab_pth))
data_pplot$Pathway = lab_pth
fplot = dotplot_top_pathways(data_pplot)
fplot
# ggsave('Factor_annotSCE.pdf', plot=fplot, width = 10, height = 10)

#--------# variance explained plot #---------#

id50 = order(gene_variances, decreasing = TRUE)[1:50]

c_sq = rowSums(fitsce$Lambda_C^2)[id50]
p_sq = rowSums(fitsce$Lambda_N^2)[id50]
gene_labs = colnames(Ys)[id50]

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
df$component = factor(df$component, levels = c("CGamma", "Psi"))

var_explained_plot = ggplot(df, aes(x = gene, y = prop, fill = component)) +
  geom_col(width = 0.9, position = position_stack(reverse = TRUE)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0001)) +
  scale_fill_manual(name = "", values = c(CGamma = "#6A3D9A", Psi = "#33A02C"),
                    labels = c(expression(C~Gamma), expression(Psi))) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 30), plot.margin = margin(8, 10, 5, 20),
        axis.text.y = element_text(size = 15))

var_explained_plot

# ggsave('var_explained_sce.pdf', plot=var_explained_plot, width = 15, height = 7.5)
