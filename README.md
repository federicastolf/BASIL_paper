# BASIL paper

This repository contains the code to replicate the simulation studies and the data
analyses in the paper *Pathway-based Bayesian factor models for omics data*
([arXiv:2601.13419](https://arxiv.org/abs/2601.13419)).

## `BASIL` package

The *R* code uses the [BASIL](https://github.com/federicastolf/BASIL) R package, which
implements a Bayesian Analysis with gene-Sets Informed Latent space (BASIL). To install
the package (in *R*):
```{r}
# if devtools is not installed yet:
# install.packages("devtools")
library(devtools)
install_github("federicastolf/BASIL")
```

## Contents

### `functs/` — shared code
- [helper.R](functs/helper.R): shared utility functions
- [FACTOR_CODE_update.R](functs/FACTOR_CODE_update.R): ROTATE
  ([Ročková & George, 2016](https://doi.org/10.1080/01621459.2015.1100620)), from
  <http://veronikarock.com/FACTOR_ANALYSIS.zip>
- [face.R](functs/face.R): FACE (covariance meta-regression with CUSP prior), from
  <https://github.com/betsybersson/covarianceMetaRegression>

### `simulations/`
- [simulations_main.R](simulations/simulations_main.R): main study — covariance and
  latent-dimension accuracy, runtime, variance ratio, and coverage (BASIL vs ROTATE vs PLIER)
- [simulations_heteroscedastic.R](simulations/simulations_heteroscedastic.R) and
  [simulations_high_noise.R](simulations/simulations_high_noise.R): same settings with
  gene-specific residual variances and with a larger residual variance, respectively
- [sim_NegBin.R](simulations/sim_NegBin.R): data generated under a negative binomial model
- [Face_sim.R](simulations/Face_sim.R): comparison against FACE
- [spectra_sim.R](simulations/spectra_sim.R) and [spectraFit_sim.py](simulations/spectraFit_sim.py): comparison against Spectra

### Data analyses
- [getGFdata.R](getGFdata.R): file to obtain the global fever data (GSE211567)
- [GFeverAnalysis.R](GFeverAnalysis.R): file to replicate the analysis on global fever data
- [scRNAseq_data.R](scRNAseq_data.R): file to replicate the scRNA-seq data analysis
- [sensitivity_WB.R](sensitivity_WB.R): sensitivity analyses on the whole-blood data 
