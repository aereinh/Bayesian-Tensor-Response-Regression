source('main_functions.R')

### Simulate data from a longitudinal model
set.seed(1234)

# Set covariates, ID/Visit/Time indices, size of tensor outcome, and true tensor rank
N_obs <- 60 # number of observations
N_subj <- 15 # number of subjects
p <- c(12,12) # tensor dimensions
D <- length(p) # number of dimensions (only allowed to be D=2 or D=3)
Visit <- rep(c(1:(N_obs/N_subj))-1,N_subj) # visit indices
Time <- Visit + rnorm(length(Visit),0,.1) # exact timepoints can be different than Visit indices, but should follow same order for a given subject
N_visit <- max(Visit)+1 # number of visits
M <- 1 # number of time-invariant covariates with time-varying effects
Ci <- array(rbinom(N_subj*M,1,.5),dim=c(N_subj,M)) # time-invariant covariates with time-varying effects
S <- 2 # number of time-invariant covariates with time-invariant effects
Xi <- array(rnorm(N_subj*S),dim=c(N_subj,S)) # covariates with time-invariant effects
K <- 1 # number of time-varying covariates with time-invariant effects
Zti <- array(rbinom(N_obs*K,1,.5),dim=c(N_obs,K)) # time-varying covariates with time-invariant effects
ID <- sort(rep(c(1:N_subj),N_visit)) # subject indices
R_true <- 2 # true rank used to generate coefficients

# Generate true model coefficients from low-rank tensor decompositions
if (D==3) {
  # generate population-level intercept using normally-distributed tensor margins
  Mu_true <- c(TP.rankR(list(array(rnorm(p[1]*R_true),dim=c(p[1],R_true)),
                             array(rnorm(p[2]*R_true),dim=c(p[2],R_true)),
                             array(rnorm(p[3]*R_true),dim=c(p[3],R_true)))))
  
  # generate subject-level intercept using normally-distributed tensor margins
  Bi_true <- t(sapply(1:N_subj,function(x) c(TP.rankR(list(array(rnorm(p[1]*R_true),dim=c(p[1],R_true)),
                                                       array(rnorm(p[2]*R_true),dim=c(p[2],R_true)),
                                                       array(rnorm(p[3]*R_true),dim=c(p[3],R_true)))))))
  
  # generate population-level time slope using normally-distributed tensor margins
  Gamma_true <- c(TP.rankR(list(array(rnorm(p[1]*R_true),dim=c(p[1],R_true)),
                                array(rnorm(p[2]*R_true),dim=c(p[2],R_true)),
                                array(rnorm(p[3]*R_true),dim=c(p[3],R_true)))))
  
  # generate subject-level time slope using normally-distributed tensor margins
  Thetai_true <- t(sapply(1:N_subj,function(x) c(TP.rankR(list(array(rnorm(p[1]*R_true),dim=c(p[1],R_true)),
                                                           array(rnorm(p[2]*R_true),dim=c(p[2],R_true)),
                                                           array(rnorm(p[3]*R_true),dim=c(p[3],R_true)))))))
  
  # generate Btm using binomially-distributed tensor margins
  Btm_true <- array(dim=c(N_visit,M,prod(p)))
  for (t in 1:N_visit) {
    for (m in 1:M) {
      Btm_true[t,m,] <- TP.rankR(
        list(
          array(rbinom(p[1]*R_true,1,.5),dim=c(p[1],R_true)),
          array(rbinom(p[2]*R_true,1,.5),dim=c(p[2],R_true)),
          array(rbinom(p[3]*R_true,1,.5),dim=c(p[3],R_true))
        )
      )
    }
  }
  
  # generate Ds using binomially-distributed tensor margins
  Ds_true <- t(sapply(1:S,function(x) c(TP.rankR(list(array(rbinom(p[1]*R_true,1,.5),dim=c(p[1],R_true)),
                                                      array(rbinom(p[2]*R_true,1,.5),dim=c(p[2],R_true)),
                                                      array(rbinom(p[3]*R_true,1,.5),dim=c(p[3],R_true)))))))
  
  # generate Gammak using binomially-distributed tensor margins
  Gammak_true <- t(sapply(1:K,function(x) c(TP.rankR(list(array(rbinom(p[1]*R_true,1,.5),dim=c(p[1],R_true)),
                                                          array(rbinom(p[2]*R_true,1,.5),dim=c(p[2],R_true)),
                                                          array(rbinom(p[3]*R_true,1,.5),dim=c(p[3],R_true)))))))
  
} else if (D==2) {
  Mu_true <- c(TP.rankR(list(array(rnorm(p[1]*R_true),dim=c(p[1],R_true)),
                             array(rnorm(p[2]*R_true),dim=c(p[2],R_true)))))
  Bi_true <- t(sapply(1:N_subj, function(x) c(TP.rankR(list(array(rnorm(p[1]*R_true),dim=c(p[1],R_true)),
                                                        array(rnorm(p[2]*R_true),dim=c(p[2],R_true)))))))
  Gamma_true <- c(TP.rankR(list(array(rnorm(p[1]*R_true),dim=c(p[1],R_true)),
                                array(rnorm(p[2]*R_true),dim=c(p[2],R_true)))))
  Thetai_true <- t(sapply(1:N_subj, function(x) c(TP.rankR(list(array(rnorm(p[1]*R_true),dim=c(p[1],R_true)),
                                                            array(rnorm(p[2]*R_true),dim=c(p[2],R_true)))))))
  Btm_true <- array(dim=c(N_visit,M,prod(p)))
  for (m in 1:M) {
    Btm_true[,m,] <-  t(sapply(1:N_visit, function(x) c(TP.rankR(list(array(rbinom(p[1]*R_true,1,.5),dim=c(p[1],R_true)),
                                                                 array(rbinom(p[2]*R_true,1,.5),dim=c(p[2],R_true)))))))
  }
  Ds_true <- t(sapply(1:S, function(x) c(TP.rankR(list(array(rbinom(p[1]*R_true,1,.5),dim=c(p[1],R_true)),
                                                       array(rbinom(p[2]*R_true,1,.5),dim=c(p[2],R_true)))))))
  Gammak_true <- t(sapply(1:K, function(x) c(TP.rankR(list(array(rbinom(p[1]*R_true,1,.5),dim=c(p[1],R_true)),
                                                           array(rbinom(p[2]*R_true,1,.5),dim=c(p[2],R_true)))))))
  
}

# Simulate first outcome Y_long from longitudinal model, with noise variance 1 across
# all observations. The `simulate_Y` function can be used to do this, which takes the tensor size,
# Time/Visit/ID indices, covariates, model coefficients, and noise variance as parameters
noise_var <- 1
Y_long <- simulate_Y(p,ID,Time,Visit,Ci,Xi,Zti,Mu=Mu_true,Bi=Bi_true,Gamma=Gamma_true,Thetai=Thetai_true,
                     Btm=Btm_true,Ds=Ds_true,Gammak=Gammak_true,noise_var)

# Fit longitudinal model on Y_long using 500 MCMC samples and rank 2. 
# Keep other parameters at the default. Note: We don't need to specify any of the null_...
# parameters, since by default the `btrr.long_rcpp` function will sample all model terms.
niter <- 500
R <- 2
btrr_results_long <- btrr.long_rcpp(Y_long,ID,Time,Visit,Ci,Xi,Zti,R,niter)

# Use the getBTRRCoef() functions to get the sampled posterior means for all model terms
Mu_est <- getBTRRCoef(btrr_results_long,term=1,burn.in=.3,find.signif=F)
plot(c(Mu_true),c(Mu_est),main="Estimated vs. True Population Intercept (M)",xlab="True",ylab="Est")
Bi_est <- getBTRRCoef(btrr_results_long,term=2,burn.in=.3,find.signif=F)
plot(c(Bi_true),c(Bi_est),main="Estimated vs. True Subject Intercept (Bi)",xlab="True",ylab="Est")
Gamma_est <- getBTRRCoef(btrr_results_long,term=3,burn.in=.3,find.signif=F)
plot(c(Gamma_true),c(Gamma_est),main="Estimated vs. True Population Time-Slope (Gamma)",xlab="True",ylab="Est")
Thetai_est <- getBTRRCoef(btrr_results_long,term=4,burn.in=.3,find.signif=F)
plot(c(Thetai_true),c(Thetai_est),main="Estimated vs. True Subject Time-Slope (Theta_i)",xlab="True",ylab="Est")

# For the terms corresponding to covariates (B_tm, D_s, and \Gamma_k), also compute the significance
# estimates using joint credible intervals (alternatively, can specify "pointwise" for standard credible intervals)
g(Btm_est,Btm_signif) %=% getBTRRCoef(btrr_results_long,term=5,burn.in=0.3,find.signif=T,signif.type="joint",alpha=0.05)
plot(c(Btm_true),c(Btm_est),col=factor(c(Btm_signif)),main="Estimated vs. True B_tm",xlab="True",ylab="Est")
g(Ds_est,Ds_signif) %=% getBTRRCoef(btrr_results_long,term=6,burn.in=0.3,find.signif=T,signif.type="joint",alpha=0.05)
plot(c(Ds_true),c(Ds_est),col=factor(c(Ds_signif)),main="Estimated vs. True D_s",xlab="True",ylab="Est")
g(Gammak_est,Gammak_signif) %=% getBTRRCoef(btrr_results_long,term=7,burn.in=0.3,find.signif=T,signif.type="joint",alpha=0.05)
plot(c(Gammak_true),c(Gammak_est),col=factor(c(Gammak_signif)),main="Estimated vs. True Gamma_k",xlab="True",ylab="Est")

# Get estimated outcome using the simulate_Y() function with estimated coefficients
Yhat_long <- simulate_Y(p,ID,Time,Visit,Ci,Xi,Zti,Mu_est,Bi_est,Gamma_est,Thetai_est,
                        Btm_est,Ds_est,Gammak_est,noise_var=0)
plot(c(Y_long),c(Yhat_long),main="Predicted vs. True Outcome",xlab="True",ylab="Predicted")

# compute MSE
MSE_long <- mean((Y_long-Yhat_long)^2)
MSE_long

# find true positive rate and true negative rate for Ds coefficient
TPR_Ds_long <- sum(Ds_signif[Ds_true!=0]==1)/sum(Ds_true!=0)
TNR_Ds_long <- sum(Ds_signif[Ds_true==0]==0)/sum(Ds_true==0)
TPR_Ds_long
TNR_Ds_long

# show traceplot for noise variance for particular observation
sigma2_mcmc <- btrr_results_long[[5]]
nMCMC <- dim(sigma2_mcmc)[2]
burn.in <- .3
obs <- 1
plot(sigma2_mcmc[obs,round(burn.in*nMCMC):nMCMC],type="l",main="Burn-in Traceplot for sigma_{ti}^2",xlab="Iteration",ylab="sigma^2")

# Calculate Deviance Information Criterion
Yvec_long <- array(Y_long,dim=c(N_obs,prod(p)))
DIC_long <- getDIC(btrr_results_long,Y_long,ID,Time,Visit,Ci,Xi,Zti,burn.in = burn.in,prog.count=10)
DIC_long$DIC1
DIC_long$DIC2











# Simulate second outcome Y_cs from cross-sectional model, with noise variance 1 across
# all observations. The `simulate_Y` function can be used to do this, with specifications
# that tell the function to exclude subject- and time-specific terms.
# Two ways of doing so are provided below.
Y_cs1 <- simulate_Y(p,ID,Time,Visit,Ci,Xi,Zti,Mu=Mu_true,Bi=0*Bi_true,Gamma=0*Gamma_true,Thetai=0*Thetai_true,
                    Btm=0*Btm_true,Ds=Ds_true,Gammak=0*Gammak_true,noise_var=0)

# repeat Xi across each subject's visits to get a N_obs x S matrix (instead of N_subj x S)
Xi_cs <- t(sapply(ID,function(x) Xi[x,]))
Y_cs2 <- simulate_Y(p,ID=NULL,Time=NULL,Visit=NULL,Ci=NULL,Xi=Xi_cs,Zti=NULL,Mu_true,Bi_true,
                    Gamma_true,Thetai_true,Btm_true,Ds_true,Gammak_true,noise_var=0)
plot(c(Y_cs1),c(Y_cs2))

noise_var <- 1
Y_cs <- simulate_Y(p,ID=NULL,Time=NULL,Visit=NULL,Ci=NULL,Xi=Xi_cs,Zti=NULL,Mu_true,Bi_true,
                   Gamma_true,Thetai_true,Btm_true,Ds_true,Gammak_true,noise_var)


# Fit cross-sectional model on Y_cs using 500 MCMC samples and rank 2. 
# Keep other parameters at the default. Note: We need to specify null_Bi, null_Gamma, 
# null_Thetai, null_Btm, and null_Gammak as FALSE to avoid sampling these (longitudinal) terms.
niter <- 500
R <- 2
btrr_results_cs <- btrr.long_rcpp(Y=Y_cs,ID=NULL,Time=NULL,Visit=NULL,Ci=NULL,Xi=Xi_cs,Zti=NULL,R=R,niter=niter,
                                  null_Bi=T, null_Gamma=T, null_Thetai=T, null_Btm=T, null_Gammak=T)

# Use the getBTRRCoef() functions to get the sampled posterior means for all model terms
Mu_est <- getBTRRCoef(btrr_results_long,term=1,burn.in=.3,find.signif=F)
plot(c(Mu_true),c(Mu_est),main="Estimated vs. True Population Intercept (M)",xlab="True",ylab="Est")

# For the terms corresponding to covariates (D_s), also compute the significance
# estimates using joint credible intervals (alternatively, can specify "pointwise" for standard credible intervals)
g(Ds_est,Ds_signif) %=% getBTRRCoef(btrr_results_long,term=6,burn.in=0.3,find.signif=T,signif.type="joint",alpha=0.05)
plot(c(Ds_true),c(Ds_est),col=factor(c(Ds_signif)),main="Estimated vs. True D_s",xlab="True",ylab="Est")

# Get estimated outcome using the simulate_Y() function with estimated coefficients
Yhat_cs <- simulate_Y(p,ID=NULL,Time=NULL,Visit=NULL,Ci=NULL,Xi=Xi_cs,Zti=NULL,
                      Mu=Mu_est,Bi=NULL,Gamma=NULL,Thetai=NULL,Btm=NULL,
                      Ds=Ds_est,Gammak=NULL,noise_var=0)
plot(c(Y_cs),c(Yhat_cs),main="Predicted vs. True Outcome",xlab="True",ylab="Predicted")

# compute MSE
MSE_long <- mean((Y_cs-Yhat_cs)^2)
MSE_long

# find true positive rate and true negative rate for Ds coefficient
TPR_Ds_long <- sum(Ds_signif[Ds_true!=0]==1)/sum(Ds_true!=0)
TNR_Ds_long <- sum(Ds_signif[Ds_true==0]==0)/sum(Ds_true==0)
TPR_Ds_long
TNR_Ds_long

# show traceplot for noise variance for particular observation
sigma2_mcmc <- btrr_results_cs[[5]]
nMCMC <- dim(sigma2_mcmc)[2]
burn.in <- .3
obs <- 1
plot(sigma2_mcmc[obs,round(burn.in*nMCMC):nMCMC],type="l",main="Burn-in Traceplot for sigma_{ti}^2",xlab="Iteration",ylab="sigma^2")

# Calculate Deviance Information Criterion
Yvec_cs <- array(Y_cs,dim=c(N_obs,prod(p)))
DIC_cs <- getDIC(btrr_results_cs,Y_cs,ID=NULL,Time=NULL,Visit=NULL,Ci=NULL,Xi=Xi_cs,Zti=NULL,burn.in=burn.in,prog.count=10)
DIC_cs$DIC1
DIC_cs$DIC2
