# BTRR
Functions to implement longitudinal or cross-sectional Bayesian Tensor Response Regression Model (BTRR) in R.

This repository contains R scripts for implementing the Gibbs sampling algorithm used to estimate model parameter posteriors in BTRR. In BTRR, each model parameter pertains to effects of scalar covariates on a tensor outcome of arbitrary dimension (most often 2 or 3), using low-rank tensor decomposition to account for spatial structure across tensor elements (e.g. effects on different voxels in a brain scan). The general regression model for the longitudinal case is as follows:
$\mathcal{Y}_{ti}=\mathcal{M}+B_i+(\Gamma+\Theta_i)\times \mathcal{T}_{ti}+\sum_{m=1}^M \mathcal{B}_{tm} c_{i,m}+\sum_{s=1}^S x_{i,s} + \sum_{k=1}^K \Gamma_k z_{ti,k} + \epsilon_{ti}$

Priors and derived conditional posteriors for each model parameter can be found in the paper "Bayesian Longitudinal Tensor Response Regression for Modeling
Neuroplasticity" (Kundu, Reinhardt, et al., 2022).

The file main_functions.R contain the R function used to implement the Gibbs sampling algorithm for l-BTRR and cs-BTRR. The function btrr.long_rcpp() repeatedly calls on btrr.1var.1visit.1iter_rcpp() across each covariate and each Gibbs sampling iteration. The resultant sampled tensor coefficients are returned, and post-burn-in samples are used to find the point estimate (mean) and 95% credible interval for each coefficient element. Functions vox.ols.reg() and vox.ols.reg.ri() are used to find the Ordinary Least Squares (OLS) estimates across each tensor element separately. Due to the inherent structure of many biological systems, BTRR may be preferable over OLS in terms of model generalizability and interpretability, for example in neuroimaging and genomics applications.

The file BTRR_demo_markdown.Rmd is a Rmarkdown file (knitted to pdf BTRR_demo_markdown.pdf) containing a brief demonstration for how to use the BTRR.R functions on simulated, longitudinal data. In this example, low-rank tensor decomposition is used to generate the true model coefficients, and resultant coefficients indicate that BTRR outperforms OLS when the low-rank generation is used. Figures showing model performance and comparison to tradiational methods (OLS) are presented here as well.
