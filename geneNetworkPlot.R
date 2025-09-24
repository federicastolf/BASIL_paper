library(PLIER)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(pheatmap)
library(scico)
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
Y = scale(Y)

source('helper_functions.R')
svd_Y = svd(Y)
plot(cumsum(svd_Y$d^2) / sum(svd_Y$d^2))
k_max = min(which((cumsum(svd_Y$d^2) / sum(svd_Y$d^2) > 0.8)))
k_max
fitBL = compute_point_estimates(Y,  allPaths, k=k_max)
Lambda_hat = fitBL$Lambda_C + fitBL$Lambda_N
CovBasil = Lambda_hat %*% t(Lambda_hat) + fitBL$sigma_sq*diag(ncol(Y))
corrBasil = cov2cor(CovBasil)

gene_variances <- apply(t(dataWholeBlood), 2, var)
id <- order(gene_variances, decreasing = TRUE)[1:100]
sub_data <- Y[,id]

basil_loadings_samples = compute_posterior_samples_cc(
  Y, fitBL$Lambda_C, fitBL$Lambda_N, fitBL$tau_C, fitBL$tau_N,
  fitBL$sigma_sq, fitBL$P_C, v0=1, sigma_sq_0=1, n_MC=500
)

cor_est = compute_correlation_posterior_samples_cc(
  basil_loadings_samples$Lambda_samples[id,,], basil_loadings_samples$sigma_sq_samples, samples=T
)

corrBasil_mean = cor_est$posterior_mean
corrBasil_qs = apply(cor_est$posterior_samples, c(1,2), function(x)(quantile(x, probs=c(0.025, 0.975))))
dim(corrBasil_qs)
corrBasil_mean[((corrBasil_qs[1,,] <0) & (corrBasil_qs[2,,] >0))] = 0
mean(corrBasil_mean==0)


row.names(corrBasil_mean) = colnames(sub_data)
colnames(corrBasil_mean) = colnames(sub_data)

library(reshape2)
library(dplyr)
df_long <- melt(corrBasil_mean[]) %>%
  filter(value < 0.9999)
df_sorted <- df_long %>%
  arrange(desc(value))

head(df_sorted, n=20)

library(tidyverse)
library(igraph)
library(ggraph)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)



rn <- corrBasil_mean[[1]]
mat <- as.matrix(corrBasil_mean)

mode(mat) <- "numeric"
C <- (mat + t(mat))/2
diag(C) <- 1
genes <- rownames(C)

A <- C
diag(A) <- 0
g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- delete_edges(g, E(g)[weight == 0])
E(g)$w_abs <- abs(E(g)$weight)
set.seed(42)
wc <- cluster_walktrap(g, weights = E(g)$w_abs)
membership <- membership(wc)
V(g)$cluster <- factor(membership)
V(g)$strength <- graph.strength(g, weights = E(g)$w_abs)
V(g)$degree <- degree(g)


dist_abs <- as.dist(1 - abs(C))
hc <- hclust(dist_abs, method = "average")
order_idx <- hc$order
C_ord <- C[order_idx, order_idx]
genes_ord <- rownames(C_ord)
cluster_vec <- membership[genes_ord]

E(g)$sign  <- ifelse(E(g)$weight > 0, "positive", "negative")
E(g)$w_abs <- abs(E(g)$weight)

g_pos <- g
E(g_pos)$weight <- E(g)$w_abs
E(g_pos)$sign   <- E(g)$sign
E(g_pos)$w_abs  <- E(g)$w_abs
V(g_pos)$cluster  <- V(g)$cluster
V(g_pos)$strength <- V(g)$strength

set.seed(123)
lay <- ggraph::create_layout(g_pos, layout = "fr")

n_labels <- min(20, igraph::gorder(g_pos))
label_genes <- names(sort(igraph::graph.strength(g_pos, weights = E(g_pos)$w_abs),
                          decreasing = TRUE))[seq_len(n_labels)]

lay$cluster  <- igraph::V(g_pos)$cluster[match(lay$name, igraph::V(g_pos)$name)]
lay$strength <- igraph::V(g_pos)$strength[match(lay$name, igraph::V(g_pos)$name)]
lay$label    <- ifelse(lay$name %in% label_genes, lay$name, NA_character_)

clu_levels <- levels(igraph::V(g_pos)$cluster)
clu_cols <- setNames(RColorBrewer::brewer.pal(max(3, length(clu_levels)) %/% 1 + 2, "Set3")[seq_along(clu_levels)], clu_levels)
edge_cols <- c(positive = "#D73027", negative = "#4575B4")


dev.new()
ggraph(lay) +
  geom_edge_link(aes(edge_alpha = w_abs, edge_width = w_abs, edge_color = sign)) +
  scale_edge_width(range = c(0.2, 0.9), guide = "none") +
  scale_edge_alpha(range = c(0.2, 0.9), guide = "none") +
  scale_edge_color_manual(values = edge_cols, name = "Edge sign") +
  geom_node_point(aes(color = cluster, size = strength)) +
  scale_color_manual(values = clu_cols, name = "Cluster") +
  scale_size_continuous(name = "Node strength", range = c(1, 3)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("Gene co-expression network")  +
  theme(legend.position = "none")
dev.off()



