
library(PLIER)
library(tidyverse)
library(reshape2)
library(gridExtra)
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

#----# fit PLIER #-----#
PLresult = PLIER(t(Y), allPaths, k = 10, doCrossval = F)
covPLIER = PLresult$Z %*% (cov(t(PLresult$B))) %*% t(PLresult$Z)
corrPLIER = cov2cor(covPLIER)

#----# fit BASIL #-----#
source('helper_functions.R')
fitBL = compute_point_estimates(Y,  allPaths, k=10)
Lambda_hat = fitBL$Lambda_C + fitBL$Lambda_N
CovBasil = Lambda_hat %*% t(Lambda_hat) + fitBL$sigma_sq*diag(ncol(Y))
corrBasil = cov2cor(CovBasil)

# random select 50 genes
set.seed(139)
id = sample(c(1:NCOL(Y)), 50, replace = F)
sub_data = Y[,id]


#------------------------# correlation Plot #---------------------#

# empirical correlation 
corP = cor(sub_data)
cor_df = melt(corP)

EmpPlot = ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue2", mid = "white", high = "red", midpoint = 0,
                       limits=c(-1,1)) +
  theme_minimal() +
  ggtitle("Empirical") +
  scale_y_discrete(limits = rev(levels(cor_df$Var1))) +
  labs(x = "", y = "", fill = "Correlation") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

EmpPlot

sub_corrBasil = corrBasil[id,id]
cor_BS = melt(sub_corrBasil)
BSPlot = ggplot(cor_BS, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  ggtitle("BASIL") +
  scale_fill_gradient2(low = "blue2", mid = "white", high = "red", midpoint = 0,
                       limits=c(-1,1)) +
  theme_minimal() +
  scale_y_discrete(limits = rev(levels(cor_BS$Var1))) +
  labs(x = "", y = "", fill = "Correlation") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

BSPlot

sub_corrPLIER = corrPLIER[id,id]
cor_PL = melt(sub_corrPLIER)

PLPlot = ggplot(cor_PL, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  ggtitle("PLIER") +
  scale_fill_gradient2(low = "blue2", mid = "white", high = "red", midpoint = 0, 
                       limits=c(-1,1)) +
  theme_minimal() +
  scale_y_discrete(limits = rev(levels(cor_PL$Var1))) +
  labs(x = "", y = "", fill = "Correlation") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")

PLPlot

CorrNP = grid.arrange(EmpPlot, BSPlot, PLPlot, nrow=1)
# ggsave('CorrWholeBLOOD.png', plot=CorrNP, device = 'png', width = 9.1, height = 3.3)


##-------------------## observed vs predicted plot ##-------------------------##


df_basil = data.frame(Observed = corP[lower.tri(corP)], 
                      Predicted = sub_corrBasil[lower.tri(sub_corrBasil)], 
                      Method = "BASIL")
df_plier = data.frame(Observed = corP[lower.tri(corP)], 
                       Predicted = sub_corrPLIER[lower.tri(sub_corrPLIER)],
                       Method = "PLIER")
lims = range(c(df_basil$Observed, df_basil$Predicted, 
               df_plier$Observed, df_plier$Predicted))

# BASIL plot
p1 = ggplot(df_basil, aes(x = Observed, y = Predicted)) +
  geom_point(color = "darkblue", shape = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(title = "BASIL", x = "Observed", y = "Predicted") +
  xlim(lims) + ylim(lims) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.title=element_text(size=16))
p1

# PLIER plot
p2 = ggplot(df_plier, aes(x = Observed, y = Predicted)) +
  geom_point(color = "darkblue", shape = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
  labs(title = "PLIER", x = "Observed", y = "Predicted") +
  xlim(lims) + ylim(lims) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=20),
        axis.title=element_text(size=16))
p2

obspPLot = grid.arrange(p1,p2, nrow=1)

# ggsave('ObsPWholeBLOOD.png', plot=obspPLot, device = 'png', width = 8, height = 4)
