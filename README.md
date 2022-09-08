# BTRR
Functions to implement longitudinal or cross-sectional Bayesian Tensor Response Regression Model (BTRR) in R.

This repository contains R scripts for implementing the Gibbs sampling algorithm used to estimate model parameter posteriors in BTRR. In BTRR, each model parameter pertains to effects of scalar covariates on a tensor outcome of arbitrary dimension (most often 2 or 3), using low-rank tensor decomposition to account for spatial structure across tensor elements (e.g. effects on different voxels in a brain scan). The general regression model for the longitudinal case is as follows:

![equation](https://latex.codecogs.com/svg.image?\mathcal{Y}_{ti}=\mathcal{M}&plus;B_i&plus;(\Gamma&plus;\Theta_i)\times&space;\mathcal{T}_{ti}&plus;\sum_{m=1}^M&space;\mathcal{B}_{tm}&space;c_{i,m}&plus;\sum_{s=1}^S&space;\mathcal{D}_s&space;x_{i,s}&plus;\sum_{k=1}^K&space;\Gamma_k&space;z_{ti,k}&plus;\epsilon_{ti})

![equation](https://latex.codecogs.com/svg.image?\epsilon_{ti}\sim&space;\mathcal{N}(0,\sigma_{ti}^2))

Priors and derived conditional posteriors for each model parameter can be found in the paper "Bayesian Longitudinal Tensor Response Regression for Modeling
Neuroplasticity" (Kundu, Reinhardt, et al., 2022).

The file main_functions.R contain the R function used to implement the Gibbs sampling algorithm for l-BTRR and cs-BTRR. The function btrr.long_rcpp() repeatedly calls on btrr.1var.1visit.1iter_rcpp() across each covariate and each Gibbs sampling iteration. The resultant sampled tensor margins (used to reconstruct the tensor coefficients) are returned, and post-burn-in samples are used to find the point estimate (mean) and credible interval for each coefficient element using the getBTRRCoef() function.

The file sample_runner.R is a R script containing a brief demonstration for how to use the BTRR functions on simulated longitudinal and cross-sectional data. In this example, low-rank tensor decomposition is used to generate the true model coefficients, and outcomes are generated from the Model using the simulate_Y() function. Figures showing the accuracy of significance estimates are presented.

**Usage (w/ list of inputs)**:

-Run MCMC using *btrr.long_rcpp(...)*, e.g. btrr_results <- btrr.long_rcpp(Y,...)
  - Y ~ Tensor-valued outcome of dimension N x p1 x p2 x p3 (3D tensor) or N x p1 x p2 (2D tensor), where N is number of observations and p=(p1,p2,p3) or (p1,p2) is the tensor (image) size. No Default.
  - ID, Time, Visit: Vectors of length N with the subject index (ID), the timepoints (Time), and the time/visit index (Visit). ID should range from 1 to number of unique subjects (N_subj), the order of values in Time and Visit for a given subject should be the same, and the minimum Visit value should be 0. Default NULL for each
  - Ci, Xi, Zti: Matrices of size N_subj x M, N_subj x S, and N x K containing covariate information corresponding to $\mathcal{B}_{tm}$, $\mathcal{D}_s$, and $\Gamma_k$, respectively. Default of NULL for each
  - R: Rank of tensor decomposition used for sampling model coefficients. Default 2
  - niter: Number of MCMC samples. Default 1000
  - null_M, null_Bi, null_Gamma, null_Thetai, null_Btm, null_Ds, null_Gammak, null_baseline: Specify whether to exclude corresponding term from model fitting. Default of FALSE for each
  - a.sig, b.sig, a.tau, b.tau, a.lam, b.lam, a.alpha, b.alpha: Hyperparameters used to define priors. Default of 1 for each
  - sigma_log_alpha: Parameter used to define proposal density variance for updating the correlation parameter $\alpha_{dr}$. Default 0.01
  - alpha.init: Initial value of $\alpha_{dr}$. Default 10
  - show.prog: Whether or not to show periodic progress updates. Default TRUE
  - prog.count: Number of samples until progress update is displayed. Default 10
  - show.allsteps: Whether or not to display updates when each model parameter is being sampled. Default FALSE
- Extract model coefficients using *getBTRRCoef(...)*
  - btrr_results ~ Outputted list from *btrr.long_rcpp()*, which contains sampled tensor margins for each of the included model terms. No Default.
  - term ~ Specifies which model term to extract the coefficient for; 1 corresponds to $\mathcal{M}$, 2 to $B_i$, 3 to $\Gamma$, 4 to $\Theta_i$, 5 to $\mathcal{B}_{tm}$, 6 to $\mathcal{D}_s$, and 7 to $\Gamma_k$. Default 6
  - burn.in ~ Specifies the proportion of starting MCMC samples to exclude from computation. Default 0.3
  - find.signif ~ Whether or not to compute significance estimates using credible intervals. Default TRUE
  - signif.type ~ Specifies the type of credible intervals which are computed, either "joint" or "pointwise". Default "joint"
  - alpha ~ Type I to be used during credible interval computation. Default 0.05
  - median ~ Whether or not to use median (vs. mean) of the sampled posterior to find the coefficient estimates. Default FALSE
  - output_mcmc ~ Whether or not to output the full post-burn-in MCMC samples for the tensor coefficient. Default FALSE
  - missing_vox ~ Vectorized indices of missing voxels which tells the function to ignore those voxels when computing point and significance estimates. Default NULL
- Use *simulate_Y()* to get predicted outcome tensors
  - p ~ Outcome tensor dimension
  - ID, Time, Visit, Ci, Xi, Zti ~ Same inputs as in *btrr.long_rcpp()*
  - Mu ~ Vectorized version of intercept $\mathcal{M}$ tensor. Default NULL
  - Bi ~ Matrix of size N_subj x prod(p) with each row containing the vectorized version of subject intercept $B_i$ for the $i$th subject. Default NULL
  - Gamma ~ Vectorized version of time-slope $\Gamma$. Default NULL
  - Thetai ~ Matrix of size N_subj x prod(p) with each row containing the vectorized version of subject time-slope $\Theta_i$ for the $i$th subject. Default NULL
  - Btm ~ Array of size N_visits x M x prod(p), where Btm[t,m,] is the vectorized version of $\mathcal{B}_{tm}$ for the $t$th visit and $m$th covariate. Default NULL
  - Ds ~ Matrix of size S x prod(p), where each row is the vectorized version of $\mathcal{D}_s$ for the $s$th covariate. Default NULL
  - Gammak ~ Matrix of size K x prod(p), where each row is the vectorized version of $\Gamma_k$ for the $k$th covariate. Default NULL
  - noise_var ~ Default 0 (keep default when finding estimated outcome; set as wish for simulated outcomes)

