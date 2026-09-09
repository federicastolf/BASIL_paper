
library(PLIER)
# library(devtools)
# install_github("federicastolf/BASIL")
library(BASIL)
library(ggplot2)
library(latex2exp)

rm(list=ls())

source("functs/helper.R")

#------------------------------------------------------------------------------#
#------------# accuracy covariance and k simulations (Fig 2a, 2c) #------------#

Nsim = 25


# Setting 1: High biological signal, p=3000
param1 = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7, 
              sd_psi = 0.1, heteroscedastic=T)
df_high_p3000_het = run_simulation_study(param1, scenario_name = "high", Nsim = Nsim,
                                         seed = 463)

# Setting 2: Low biological signal, p=3000
param2 = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.4, 
              sd_psi = 0.7, heteroscedastic=T)
df_low_p3000_het = run_simulation_study(param2, scenario_name = "low", Nsim = Nsim, 
                                        seed = 463)

# Setting 3: High biological signal, p=1000
param3 = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7, 
              sd_psi = 0.1, heteroscedastic=T)
df_high_p1000_het = run_simulation_study(param3, scenario_name = "high", Nsim = Nsim,
                                         seed = 463)

# Setting 4: Low biological signal, p=1000
param4 = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.4, 
              sd_psi = 0.7, heteroscedastic=T)
df_low_p1000_het = run_simulation_study(param4,scenario_name = "low", Nsim = Nsim, 
                                        seed = 463)


# all results
Simboxplot_df_het = rbind(df_high_p3000_het, df_low_p3000_het, df_high_p1000_het, df_low_p1000_het)


#------# MSE boxplot #-------#

nl = c("1000"="1000 genes", "3000"="3000 genes")

# N.B. filter out BASIL_posterior or decide if you want to keep it 

Fnplot = ggplot(Simboxplot_df_het |>
                  dplyr::filter(model %in%  c('PLIER', 'BASIL_posterior', 'ROTATE')), aes(x = scenario, y = err_norm, fill = model))+
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(
    values = c(
      "PLIER"   = "#c85200", 
      "BASIL_posterior" = "#009E73",
      "ROTATE"          = "#1170aa"
    ),
    breaks = c("BASIL_posterior", "PLIER", "ROTATE"),
    labels = c("BASIL", "PLIER", "ROTATE")
  ) +
  facet_wrap(~ p, scales = "fixed", labeller = as_labeller(nl)) +
  xlab("Biological signal") + ylab("Error") +
  theme_light() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=15),legend.key.size = unit(1,"line"),
        legend.box.spacing = unit(0.1,"line"),
        strip.text = element_text(size = 16, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=15),)
Fnplot

# ggsave(filename = "results_sim/Errplot_het.png", plot=Fnplot,  width = 9, height = 5)

latent_factors_plot = ggplot(Simboxplot_df_het |>
                               dplyr::filter(model %in%  c('PLIER-true k', 'BASIL_posterior', 'ROTATE')), aes(x = scenario, y = err_factors, fill = model))+
  geom_boxplot(alpha=0.7) +
  # scale_fill_manual(values = c("green3","red", "steelblue")) +
  #scale_fill_manual(values =c("#009E73", "#c85200","#1170aa")) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "PLIER-true k"   = "#c85200", 
      "BASIL_posterior" = "#009E73",
      "ROTATE"          = "#1170aa"
    ),
    breaks = c("BASIL_posterior", "PLIER-true k", "ROTATE"),
    labels = c("BASIL", "PLIER", "ROTATE")
  ) +
  facet_wrap(~ p, scales = "fixed", labeller = as_labeller(nl)) +
  xlab("Biological signal") + ylab("Error") +
  theme_light() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=15),legend.key.size = unit(1,"line"),
        legend.box.spacing = unit(0.1,"line"),
        strip.text = element_text(size = 16, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=15),)
latent_factors_plot

# ggsave(filename = "results_sim/ErrLatentFactorsplot_het.png", plot=latent_factors_plot,  width = 9, height = 5)



#------------------------------------------------------------------------------#
#-----------------# Uncertainty quantification (Fig 2d) #----------------------#

subsample_size = 200
subsample_index = 1:subsample_size

# Setting 1: High biological signal, p=3000
coverage_high_p3000_het <- run_coverage_simulation(param1, scenario_name = "high", 
                                               subsample_index, alpha = 0.05, Nsim = Nsim, seed = 463)

# Setting 2: Low biological signal, p=3000
coverage_low_p3000_het <- run_coverage_simulation(param2, scenario_name = "low", 
                                              subsample_index, alpha = 0.05, Nsim = Nsim, seed = 463)

# Setting 3: High biological signal, p=1000
coverage_high_p1000_het <- run_coverage_simulation(param3, scenario_name = "high", 
                                               subsample_index, alpha = 0.05, Nsim = Nsim, seed = 463)

# Setting 4: Low biological signal, p=1000
coverage_low_p1000_het <- run_coverage_simulation(param4, scenario_name = "low", 
                                              subsample_index, alpha = 0.05, Nsim = Nsim, seed = 463)

SimUQ_het = rbind(coverage_high_p3000_het, coverage_low_p3000_het, coverage_high_p1000_het,
              coverage_low_p1000_het)

UQplot_het = ggplot(SimUQ_het, aes(x = scenario, y = coverage))+
  geom_boxplot(alpha=0.7, fill="lightblue") +
  facet_wrap(~ p, scales = "fixed", labeller = as_labeller(nl)) +
  geom_hline(aes(yintercept=0.95), linetype = "dashed") +
  xlab("Biological signal") +
  ylab("Coverage") +
  theme_light() +
  ylim(c(0.6,1)) +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=15),legend.key.size = unit(1,"line"),
        legend.box.spacing = unit(0.1,"line"),
        strip.text = element_text(size = 16, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=15),)
UQplot_het
# ggsave(filename = "results_sim/Coverageplot_het.png", plot=UQplot_het,  width = 5, height = 5)




Simboxplot_df_het |>
  dplyr::filter(model %in%  c('PLIER', 'BASIL', 'ROTATE'))  |>
  dplyr::group_by(model, scenario, p) |>
  dplyr::summarise(
    mean = mean(k_est, na.rm = TRUE),
    sd   = sd(k_est, na.rm = TRUE),
    .groups = "drop"
  )



