# BASIl-paper

This repository contains the R code to replicate the simulations study and the data analysis in the paper Pathway-based Bayesian factor models for gene expression data.

## `BASIL` package

The *R* code uses the [BASIL](https://github.com/federicastolf/BASIL) R package which implements a Bayesian Analysis with gene-Sets Informed Latent space (BASIL). To install the package (in *R*):
```{r}
# if devtools is not installed yet:
# install.packages("devtools")
library(devtools)
install_github("federicastolf/BASIL")
```
## Contents

The repository includes the following files:
- [simulations](https://github.com/federicastolf/BASIL_paper/blob/main/simulations.R) file to replicate the simulation study
- [getGFdata](https://github.com/federicastolf/BASIL_paper/blob/main/getGFdata.R) file to obtain the global fever data (GSE211567)
- [GFanalysis](https://github.com/federicastolf/BASIL_paper/blob/main/GFeverAnalysis.R) file to replicate the analysis on global fever data
