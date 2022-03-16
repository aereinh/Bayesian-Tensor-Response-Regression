# BTRR
Functions to implement Bayesian Tensor Response Regression Model (BTRR) in R

This repository contains R scripts for implementing the Gibbs sampling algorithm used to estimate model parameter posteriors in BTRR. In BTRR, each model parameter pertains to effects of scalar covariates on a tensor outcome of arbitrary dimension (most often 2 or 3), using low-rank tensor decomposition to account for spatial structure across tensor elements (e.g. effects on different voxels in a brain scan). The overall regression model, along with the priors and conditional posteriors for each model parameter, are found in Posterior_derivations.pdf.

The file BTRR.R contains R functions used to implement the Gibbs sampling algorithm. The function tensor.reg() repeatedly calls on tensor.reg.1var.1iter() across each covariate and each MCMC iteration. The resultant sampled coefficients are returned, and post-burn-in samples are used to find the point estimate (mean) and 95% credible interval for each coefficient element. Functions vox.ols.reg() and vox.ols.reg.ri() are used to find the Ordinary Least Squares (OLS) estimates across each tensor element separately. Due to the inherent structure of many biological systems, BTRR may be preferable over OLS in terms of model generalizability and interpretability.

The file demo.R contains a brief demonstration for how to use the BTRR.R functions on simulated, longitudinal data. In this example, low-rank tensor decomposition is used to generate the true model coefficients, and resultant coefficients indicate that BTRR outperforms OLS when the low-rank generation is used.
