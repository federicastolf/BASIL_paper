library(PLIER)

rm(list=ls())

source('helper_functions.R')
source("simulation_helper.R")


#------------------------------------------------------------------------------#
#------------# accuracy covariance and k simulations (Fig 2a, 2c) #------------#

Nsim = 25

# Setting 1: High biological signal, p=3000
param1 = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7, 
              sd_psi = 0.1)
df_high_p3000 = run_simulation_study(param1, scenario_name = "high", Nsim = Nsim,
                                     seed = 463)

# Setting 2: Low biological signal, p=3000
param2 = list(n = 500, p = 3000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.4, 
              sd_psi = 0.7)
df_low_p3000 = run_simulation_study(param2, scenario_name = "low", Nsim = Nsim, 
                                     seed = 463)

# Setting 3: High biological signal, p=1000
param3 = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7, 
              sd_psi = 0.1)
df_high_p1000 = run_simulation_study(param3, scenario_name = "high", Nsim = Nsim,
                                      seed = 463)

# Setting 4: Low biological signal, p=1000
param4 = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.4, 
              sd_psi = 0.7)
df_low_p1000 = run_simulation_study(param4,scenario_name = "low", Nsim = Nsim, 
                                     seed = 463)

# all results
Simboxplot_df = rbind(df_high_p3000, df_low_p3000, df_high_p1000, df_low_p1000)

#------# MSE boxplot #-------#

nl = c("1000"="1000 genes", "3000"="3000 genes")

Fnplot = ggplot(Simboxplot_df, aes(x = scenario, y = err_norm, fill = model))+
  geom_boxplot(alpha=0.7) +
  # scale_fill_manual(values = c("green3","red", "steelblue")) +
  scale_fill_manual(values =c("#009E73", "#c85200","#1170aa")) +
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

# ggsave(filename = "results_sim/Errplot.png", plot=Fnplot,  width = 9, height = 5)


#------------------------------------------------------------------------------#
#------------------------# computation time (Fig 2b) #-------------------------#

param5 = list(n = 500, p = 4000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7, 
              sd_psi = 0.1)
df_high_p4000 = run_simulation_study(param5,scenario_name = "high", Nsim = Nsim, 
                                     seed = 463)
param6 = list(n = 500, p = 2000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7, 
              sd_psi = 0.1)
df_high_p2000 = run_simulation_study(param6,scenario_name = "high", Nsim = Nsim, 
                                     seed = 463)
param7 = list(n = 500, p = 5000, k = 10, q = 500, sigma_sq_0 = 15, sd_gamma = 0.7, 
              sd_psi = 0.1)
df_high_p5000 = run_simulation_study(param7,scenario_name = "high", Nsim = Nsim, 
                                     seed = 463)

Time_temp = Simboxplot_df %>% filter(scenario=="high")
TimeSim = rbind(Time_temp, df_high_p4000, df_high_p2000, df_high_p5000)
TimeSim = TimeSim %>% dplyr::select(time, model, p)

#--# plot
Timeavg = TimeSim %>% group_by(model,p) %>% summarise(mean = mean(time))
Timeavg$p = as.numeric(Timeavg$p)
Timeplot = ggplot(Timeavg, aes(x = p, y = log(mean), color = model, group = model)) + 
  geom_point() + geom_line() +
  # scale_color_manual(values = c("green3","red", "steelblue")) +
  scale_color_manual(values = c("#009E73", "#c85200","#1170aa")) +
  xlab("Number of genes") + ylab("runtime (log scale)") +
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

Timeplot

# ggsave(filename = "results_sim/Timeplot.png", plot=Timeplot,  width = 7.5, height = 5)


#------------------------------------------------------------------------------#
#-------------# Ratio variances and biological signal (Fig 2e) #---------------#

sd_gammaL = c(0.001, 0.09, 0.14, 0.2, 0.25 , 0.28, 0.31, 0.34, 0.37)

set.seed(463)
seeds_g = sample.int(9000, Nsim)
tauC = tauN = ratio = matrix(0, length(sd_gammaL), Nsim)

for(g in 1:length(sd_gammaL)){
  param = list(n = 500, p = 1000, k = 10, q = 500, sigma_sq_0 = 2, sd_gamma = sd_gammaL[g],
               sd_psi = 1)
  
  for(s in 1:Nsim){
    #simulate data
    datas = syntheticData(n = param$n, p = param$p, k = param$k, q = param$q, 
                          sigma_sq_0 = param$sigma_sq_0, sd_gamma = param$sd_gamma,
                          sd_psi = param$sd_psi, mseed = seeds_g[s])
    Ys = datas$Y
    Cs = datas$C
    
    # compute BASIL
    fitBASIL = compute_point_estimates(Ys, Cs, k = param$k)
    tauN[g,s] = fitBASIL$tau_N 
    tauC[g,s] = fitBASIL$tau_C
    ratio[g,s] = tauC[g,s]/tauN[g,s]
  }
}

## plot
vR = c(t(ratio))
sig = rep(1:9, each = 25)
dataPlot_bs = cbind.data.frame(vR,sig)

p_bs = ggplot(dataPlot:bs, aes(y = vR, group = sig))+
  geom_boxplot(alpha=0.7, fill="lightblue1") +
  xlab("Biological signal") +
  ylab(TeX("variance $\\Gamma$ / variance $\\Psi$")) +
  geom_hline(yintercept=1, linetype="dashed") +
  theme_light() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=15),legend.key.size = unit(1,"line"),
        legend.box.spacing = unit(0.1,"line"),
        strip.text = element_text(size = 16, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=13),
        axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=16),)
p_bs

# ggsave(filename = "results_sim/LowHigh.png", plot=p_bs,  width = 8.5, height = 5)


#------------------------------------------------------------------------------#
#-----------------# Uncertainty quantification (Fig 2d) #----------------------#

subsample_index = 1:subsample_size

# Setting 1: High biological signal, p=3000
coverage_high_p3000 <- run_coverage_simulation(
  param1, subsample_index, scenario_name = "high", subsample_size = 200, alpha = 0.05, 
  Nsim = Nsim, seed = 463)

# Setting 2: Low biological signal, p=3000
coverage_low_p3000 <- run_coverage_simulation(
  param2, subsample_index, scenario_name = "low", subsample_size = 200, alpha = 0.05, 
  Nsim = Nsim, seed = 463)

# Setting 3: High biological signal, p=1000
coverage_high_p1000 <- run_coverage_simulation(
  param3, subsample_index, scenario_name = "high", subsample_size = 200, alpha = 0.05, 
  Nsim = Nsim, seed = 463)

# Setting 4: Low biological signal, p=1000
coverage_low_p1000 <- run_coverage_simulation(
  param4, subsample_index, scenario_name = "low", subsample_size = 200, alpha = 0.05, 
  Nsim = Nsim, seed = 463)

SimUQ = rbind(coverage_high_p3000$data, coverage_low_p3000$data, coverage_high_p5000$data,
  coverage_low_p5000$data)

UQplot = ggplot(SimUQ, aes(x = scenario, y = coverage))+
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
UQplot

# ggsave(filename = "results_sim/Coverageplot.png", plot=UQplot,  width = 5, height = 5)

