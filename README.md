# BTRR
Functions to implement longitudinal or cross-sectional Bayesian Tensor Response Regression Model (BTRR) in R.

This repository contains R scripts for implementing the Gibbs sampling algorithm used to estimate model parameter posteriors in BTRR. In BTRR, each model parameter pertains to effects of scalar covariates on a tensor outcome of arbitrary dimension (most often 2 or 3), using low-rank tensor decomposition to account for spatial structure across tensor elements (e.g. effects on different voxels in a brain scan). The general regression model for the longitudinal case is as follows:

![equation](https://latex.codecogs.com/svg.image?\mathcal{Y}_{ti}=\mathcal{M}&plus;B_i&plus;(\Gamma&plus;\Theta_i)\times&space;\mathcal{T}_{ti}&plus;\sum_{m=1}^M&space;\mathcal{B}_{tm}&space;c_{i,m}&plus;\sum_{s=1}^S&space;\mathcal{D}_s&space;x_{i,s}&plus;\sum_{k=1}^K&space;\Gamma_k&space;z_{ti,k}&plus;\epsilon_{ti})

![equation](https://latex.codecogs.com/svg.image?\epsilon_{ti}\sim&space;\mathcal{N}(0,\sigma_{ti}^2))

Priors and derived conditional posteriors for each model parameter can be found in the paper "Bayesian Longitudinal Tensor Response Regression for Modeling
Neuroplasticity" (Kundu, Reinhardt, et al., 2022).

The file main_functions.R contain the R function used to implement the Gibbs sampling algorithm for l-BTRR and cs-BTRR. The function btrr.long_rcpp() repeatedly calls on btrr.1var.1visit.1iter_rcpp() across each covariate and each Gibbs sampling iteration. The resultant sampled tensor margins (used to reconstruct the tensor coefficients) are returned, and post-burn-in samples are used to find the point estimate (mean) and credible interval for each coefficient element using the getBTRRCoef() function.

The file sample_runner.Rmd is a Rmarkdown file containing a brief demonstration for how to use the BTRR functions on simulated longitudinal data. In this example, low-rank tensor decomposition is used to generate the true model coefficients, and outcomes are generated from the Model using the simulate_Y() function. Figures showing the accuracy of significance estimates are presented.
