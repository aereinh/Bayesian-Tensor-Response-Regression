source('rcppFunctions.R')
source('helper_functions.R')
library(dplyr)
library(tidyr)
library(gsubfn)
library(MASS)
library(GIGrvg)
library(rlist)
library(Rcpp)

### (MAIN FUNCTION) btrr.long_rcpp: Performs MCMC sampling for all terms of longitudinal model
#
# USAGE: btrr_results <- btrr.long_rcpp(Y, ID, Time, Visit, Ci, Xi, Zti, R, niter,
#                                       null_M, null_Bi, null_Gamma, null_Thetai, null_Btm, null_Ds, null_Gammak,null_baseline,
#                                       a.sig, b.sig, a.tau, b.tau, a.lam, b.lam, a.alpha, b.alpha,sigma_log_alpha,alpha.init, 
#                                       show.prog=T, prog.count=10, show.allsteps=F)
#
#        g(beta.results,w.results,lambda.results,tau.results,sigma2.results,alpha.results) %=% btrr_results
#
# INPUTS: 
#          Y         ----> Tensor-valued outcome of dimension N x p1 x p2 x p3 (3D tensor) or N x p1 x p2 (2D tensor), 
#                          where N is number of observations and p=(p1,p2,p3) or (p1,p2) is the tensor (image) size.
#
#          ID        ----> Integer vector of length N indexing subject. ID should range from 1 to Ni.
#                          Set to NULL when modeling cross-sectional data.
# 
#          Time      ----> Numeric vector of length N representing times of each visit, where nth entry pertains to nth slice of Y.
#                          Set to NULL when modeling cross-sectional data.
# 
#          Visit     ----> Integer vector of length N indexing visits. The first/baseline visit should be coded as 0, and subsequent visit indices as 1, 2, etc.
#                          Set to NULL when modeling cross-sectional data.
# 
#          Ci        ----> Numeric matrix of dimension Ni x M, where Ni is the number of unique subjects in the dataset (1<=Ni<=N). 
#                          The ith row of Ci pertains to ith subject, and each of the M columns are covariates assumed to have a time-varying effect on Y
#                          (i.e. each of the M covariates has a different effect on each Visit). Set to NULL when modeling cross-sectional data.
# 
#          Xi        ----> Numeric matrix of dimension Ni x S, where Ni is the number of unique subjects in the dataset (1<=Ni<=N).
#                          The ith row of Xi pertains to the ith subject, and each of the S columns are covariates assumed to have a time-invariant effect on Y
#                          (i.e. the sth covariate has the same effect on Y across all Visits).
# 
#          Zti       ----> Numeric matrix of dimension N x K, where N is number of observations in the dataset. The nth row pertains to
#                          the nth slice of Y (i.e. Y[n,,,] or Y[n,,]), and each of the K columns are time-varying covariates.
#                          Set to NULL when modeling cross-sectional data.
# 
#          R         ----> Integer vector of length 7, or scalar integer specifying rank used for low-rank decomposition.
#                          If R is a vector, each element corresponds to rank of corresponding coefficient in longitudinal model term:
#                          i.e. M, B_i, Gamma, Theta_i, B_tm, D_s, and Gamma_k. If R is a scalar, all coefficients are assumed to have equal rank.
#          
#          niter     ----> Integer specofying number of MCMC samples
# 
# null_M-null_Gammak ----> Booleans specifying whether or not to exclude corresponding term from model fitting. Model terms are:
#                          M (population intercept), B_i (subject intercept), Gamma (population time-slope), Theta_i (subject time-slope),
#                          B_tm (time-varying effects), D_s (time-invariant effects), Gamma_k (effects for time-varying covariates).
#                          Cross-sectional modeling will be performed if null_Bi=null_Gamma=null_Thetai=null_Btm=null_Gammak=TRUE.
#                          
# null_baseline      ----> In longitudinal case where B_tm term is included, set null_baseline=TRUE in order to exclude the effect
#                          of Ci on the baseline visit (Visit=0). Default is FALSE.
#
# a.sig-b.alpha      ----> Numeric scalars specifying hyperparameters of model priors, namely for terms \sigma^2, \tau, \lambda, and \alpha parameters
#
# sigma_log_alpha    ----> Numeric scalar specifying the variance of proposal/jumping density for \alpha parameter, which is assumed to be Log-Normal.
#                          Too high choices of sigma_log_alpha may lead to sporadic coefficients estimates across MCMC samples
#
# alpha.init         ----> Numeric scalar specifying the initial value of \alpha when performing MCMC sampling
#
# show.prog/allsteps ----> Booleans specifying whether to periodically show which MCMC step and/or which term is being sampled
#
# prog.count         ----> Integer which specifies how often to show progress updates (i.e. how many steps until update)
#
#
#
# OUTPUTS:
#         * NOTE: All but one outputted variable contains different values for each of the 7 distinct sets of model terms.
#                Those terms are as follows:
#                   1. M ~ population-level intercept
#                   2. B_i ~ subject intercept
#                   3. \Gamma ~ population time-slope
#                   4. \Theta_i ~ subject time-slope
#                   5. B_{tm} ~ time-varying coefficients corresponding to time-invariant covariates (C_i)
#                   6. D_s ~ time-invariant coefficients corresponding to time-invariant covariates (X_i)
#                   7. Gamma_k ~ time-invariant coefficients corresponding to time-varying covariates (Z_{ti})
#                 If the inputs specify that certain parameters should be excluded from model fitting, the outputted parameters
#                 corresponding to the excluded terms will contain NAs.
#                 Denote nMCMC= Number of MCMC iterations applied + 1 initialization iteration
#
#         beta.store      ----> (List of length 7) Each element of the list contains a list of length D (i.e. # of tensor dimensions, D=2 or D=3),
#                               where each element of those sublist contains a (p[d] x R x J x nMCMC) array with the sampled tensor margins (p[d] elements) 
#                               across all ranks (R), covariates/subjects (J), and MCMC iterations (nMCMC) for that particular term
#                         
#         w.store         -----> (List of length 7) Each element of the list contains a (D x R x nMCMC) array with the diagonal entries of the margin covariance matrix (w_{dr})
#                                for the corresponding term across all MCMC iterations, where D is the dimension size and R is the applied tensor decomposition rank
#
#         lambda.store    -----> (List of length 7) Each element of the list contains a (D x R x nMCMC) array with each term's rate parameter \lambda_{dr} across all MCMC samples
#
#         tau.store       -----> (List of length 7) Each element of the list contains a matrix of size (J x nMCMC), where J is the number of subjects/covariates for the
#                                corresponding model term (J=1 for M and \Gamma). These values provide the sampled global scaling factors \tau_j
#
#         sigma2.store    -----> N x nMCMC array containing the samples of noise variance \sigma_{ti} across each observation (total of N observations)
#
#         alpha.store     -----> (List of length 7) Each element of the list contains a (D x R x nMCMC) array with each term's margin covariance 
#                                correlation parameter \alpha_{dr} across all MCMC samples
btrr.long_rcpp <- function(Y,ID=NULL,Time=NULL,Visit=NULL,Ci=NULL,Xi=NULL,Zti=NULL,R=2,niter=1000,
                           null_M=F, null_Bi=F,null_Gamma=F, null_Thetai=F, null_Btm=F, null_Ds=F, null_Gammak=F, null_baseline=F,  
                           a.sig=1, b.sig=1, a.tau=1, b.tau=1,a.lam=1, b.lam=1, a.alpha=1, b.alpha=1,sigma_log_alpha=0.01,alpha.init=10, 
                           show.prog=T, prog.count=10, show.allsteps=F) {
    nterms <- 7 # {M, B_i, \Gamma, \Theta_i, B_{tm}, D_s, \Gamma_k}
    
    # define size parameters and reformat inputs if necessary
    if (length(R)==1) {R <- rep(R,nterms)}
    N <- dim(Y)[1]
    p <- dim(Y)[-1]
    D <- length(p)
    Yvec <- arrayC(Y,dim=c(N,prod(p)))
    if (is.null(Visit)) {
      Visit <- rep(0,N)
      null_Gamma=null_Thetai=null_Btm=T
    }
    if (is.null(Time)) {
      Time <- Visit
    }
    if (is.null(ID)) {
      ID <- rep(1,N)
      null_Bi=null_Thetai=T
    }
    Ni <- length(unique(ID)); if (Ni==1) {null_Bi=null_Thetai=T}
    Nt <- length(unique(Visit)); if (Nt==1) {null_Btm=T}
    if (is.null(Ci)) {
      Ci <- as.matrix(rep(0,Ni))
      null_Btm=T
    } else if (is.vector(Ci)) {
      Ci <- as.matrix(Ci)
    }
    M <- ncol(Ci)
    if (is.null(Xi)) {
      Xi <- as.matrix(rep(0,Ni))
      null_Ds=T
    } else if (is.vector(Xi)) {
      Xi <- as.matrix(Xi)
    }
    S <- ncol(Xi)
    if (is.null(Zti)) {
      Zti <- as.matrix(rep(0,N))
      null_Gammak=T
    } else if (is.vector(Zti)) {
      Zti <- as.matrix(Zti)
    }
    K <- ncol(Zti)
    
    # initialize storage variables
    beta.store <- vector(mode="list",nterms)
    w.store <- vector(mode="list",nterms)
    lambda.store <- vector(mode="list",nterms)
    alpha.store <- vector(mode="list",nterms)
    for (ii in 1:nterms){
      beta.store[[ii]] <- vector(mode="list",D)
      w.store[[ii]] <- array(dim=c(D,R[ii],niter+1))
      w.store[[ii]][,,1] <- 1
      lambda.store[[ii]] <- array(dim=c(D,R[ii],niter+1))
      lambda.store[[ii]][,,1] <- 1
      alpha.store[[ii]] <- array(dim=c(D,R[ii],niter+1))
      alpha.store[[ii]][,,1] <- alpha.init
    }
    tau.store <- list(rep(NA,niter+1), array(dim=c(Ni,niter+1)), rep(NA,niter+1), 
                      array(dim=c(Ni,niter+1)), array(dim=c(Nt,M,niter+1)), 
                      array(dim=c(S,niter+1)), array(dim=c(K,niter+1)))
    tau.store[[1]][1]<-1; tau.store[[2]][,1]<-1; tau.store[[3]][1]<-1
    tau.store[[4]][,1]<-1; tau.store[[5]][,,1]<-1; tau.store[[6]][,1]<-1
    tau.store[[7]][,1]<-1
    
    if (null_baseline) {tau.store[[5]][1,,] <- 1}
    sigma2.store <- array(dim=c(N,niter+1)); sigma2.store[,1]<-1
    for (d in 1:D) {
      beta.store[[1]][[d]] <- array(dim=c(p[d],R[1],niter+1)) # M (intercept)
      beta.store[[1]][[d]][,,1] <- rnorm(p[d]*R[1]) # M (intercept)
      beta.store[[2]][[d]] <- array(dim=c(p[d],R[2],Ni,niter+1)) # Bi
      beta.store[[2]][[d]][,,,1] <- rnorm(p[d]*R[2]*Ni) # Bi
      beta.store[[3]][[d]] <- array(dim=c(p[d],R[3],niter+1)) # Gamma
      beta.store[[3]][[d]][,,1] <- rnorm(p[d]*R[3]) # Gamma
      beta.store[[4]][[d]] <- array(dim=c(p[d],R[4],Ni,niter+1)) # Theta_i
      beta.store[[4]][[d]][,,,1] <- rnorm(p[d]*R[4]*Ni) # Theta_i
      beta.store[[5]][[d]] <- array(dim=c(p[d],R[5],Nt,M,niter+1)) # B_tm
      beta.store[[5]][[d]][,,,,1] <- rnorm(p[d]*R[5]*Nt*M) # B_tm
      if (null_baseline) {beta.store[[5]][[d]][,,1,,] <- 0}
      beta.store[[6]][[d]] <- array(dim=c(p[d],R[6],S,niter+1)) # D_s
      beta.store[[6]][[d]][,,,1] <- rnorm(p[d]*R[6]*S) # D_s
      beta.store[[7]][[d]] <- array(dim=c(p[d],R[7],K,niter+1)) # Gamma_k
      beta.store[[7]][[d]][,,,1] <- rnorm(p[d]*R[7]*K) # Gamma_k
    }

    if (!null_M) {
      M_iter <- c(TP.rankR(lapply(beta.store[[1]],function(x) x[,,1])))
    } else {
      M_iter <- rep(0,prod(p))
    }
    if (!null_Bi) {
      Bi_iter <- t(sapply(1:Ni,function(i) TP.rankR(lapply(beta.store[[2]],function(x) x[,,i,1]))))
    } else {
      Bi_iter <- array(0,dim=c(Ni,prod(p)))
    }
    if (!null_Gamma) {
      Gamma_iter <- c(TP.rankR(lapply(beta.store[[3]],function(x) x[,,1])))
    } else {
      Gamma_iter <- rep(0,prod(p))
    }
    if (!null_Thetai) {
      Thetai_iter <- t(sapply(1:Ni,function(i) TP.rankR(lapply(beta.store[[4]],function(x) x[,,i,1]))))
    } else {
      Thetai_iter <- array(0,dim=c(Ni,prod(p)))
    }
    Btm_iter <- array(0,dim=c(Nt,M,prod(p)))
    if (!null_Btm) {
      for (tt in 1:Nt) {
        for (mm in 1:M) {
          Btm_iter[tt,mm,] <- c(TP.rankR(lapply(beta.store[[5]],function(x) x[,,tt,mm,1])))
        }
      }
    }
    Btm_iter_cond <- arrayC(Btm_iter,dim=c(Nt*M,prod(p)))
    if (!null_Ds) {
      Ds_iter <- t(sapply(1:S, function(s) TP.rankR(lapply(beta.store[[6]], function(x) x[,,s,1]))))
    } else {
      Ds_iter <- array(0,dim=c(S,prod(p)))
    }
    if (!null_Gammak) {
      Gammak_iter <- t(sapply(1:K,function(k) TP.rankR(lapply(beta.store[[7]], function(x) x[,,k,1]))))
    } else {
      Gammak_iter <- array(0,dim=c(K,prod(p)))
    }
    

    for (iter in 1:niter) {
      
      if (show.prog && iter%%prog.count==0) {
        print(iter)
      }
      
      # 1. update M terms
      if (!null_M) {
        if (show.prog && iter%%prog.count==0 && show.allsteps) {
          print(paste0(iter,': Updating M'))
        }
        Y_M_vec <- getResid_fullModel(Yvec,ID,Time,Visit,Ci,Xi,Zti,0*M_iter,Bi_iter,
                                      Gamma_iter,Thetai_iter,Btm_iter_cond,Ds_iter,Gammak_iter,M)
        
        Y_M <- arrayC(Y_M_vec, dim=c(N,p))
        g(beta.new,w.new,lambda.new,tau.new,sigma2.new,alpha.new) %=% btrr.1var.1visit.1iter_rcpp(Y = Y_M, xc = as.matrix(rep(1,N)), R=R[1], 
                                                                                                 beta.prev = lapply(beta.store[[1]],function(x) arrayC(x[,,iter],dim=c(dim(x)[1],R[1]))), w.prev=arrayC(w.store[[1]][,,iter],dim=c(D,R[1])), 
                                                                                                 lambda.prev = arrayC(lambda.store[[1]][,,iter],dim=c(D,R[1])), tau.prev = tau.store[[1]][iter], 
                                                                                                 sigma2.prev = sigma2.store[,iter], alpha.prev = arrayC(alpha.store[[1]][,,iter],dim=c(D,R[1])),
                                                                                                 update.w = T, update.lambda = T, update.tau = T, update.sigma2 = F, update.alpha = T, sigma_log_alpha = sigma_log_alpha)
        

        for (d in 1:D) {beta.store[[1]][[d]][,,iter+1] <- beta.new[[d]]}
        w.store[[1]][,,iter+1] <- w.new
        lambda.store[[1]][,,iter+1] <- lambda.new
        tau.store[[1]][iter+1] <- tau.new
        alpha.store[[1]][,,iter+1] <- alpha.new
        
        M_iter <- c(TP.rankR(beta.new))
      }
      
      # 2. update Bi terms
      if (!null_Bi) {
        if (show.prog && iter%%prog.count==0 && show.allsteps) {
          print(paste0(iter,': Updating Bi'))
        }
        c.dr_array_Bi <- array(0,dim=c(D,R[2]))
        for (i in 1:Ni){
          Visit_i <- Visit[ID==i]
          Time_i <- Time[ID==i]
          Zti_i <- as.matrix(Zti[ID==i,])
          ID_i <- ID[ID==i]
          Y_Bi_vec <- getResid_fullModel(Yvec[ID==i,],ID_i,Time_i,Visit_i,Ci,Xi,Zti_i,
                                         M_iter,0*Bi_iter,Gamma_iter,Thetai_iter,Btm_iter_cond,Ds_iter,Gammak_iter,M)
          Y_Bi <- arrayC(Y_Bi_vec, dim=c(dim(Y_Bi_vec)[1],p))
          g(beta.new, w.new, lambda.new, tau.new, sigma2.new, alpha.new) %=% btrr.1var.1visit.1iter_rcpp(Y = Y_Bi, xc = as.matrix(rep(1,dim(Y_Bi)[1])), R=R[2],
                                                                                                            beta.prev = lapply(beta.store[[2]],function(x) arrayC(x[,,i,iter],dim=c(dim(x)[1],R[2]))), w.prev=arrayC(w.store[[2]][,,iter],dim=c(D,R[2])),
                                                                                                            lambda.prev = arrayC(lambda.store[[2]][,,iter],dim=c(D,R[2])), tau.prev = tau.store[[2]][i,iter],
                                                                                                            sigma2.prev = sigma2.store[ID==i,iter], alpha.prev = arrayC(alpha.store[[2]][,,iter],dim=c(D,R[2])),
                                                                                                            update.w = F, update.lambda = F, update.tau = T, update.sigma2 = F, update.alpha = F)
          for (d in 1:D) {beta.store[[2]][[d]][,,i,iter+1] <- beta.new[[d]]}
          tau.store[[2]][i,iter+1] <- tau.new

          Bi_iter[i,] <- c(TP.rankR(beta.new))
          
          if (D==3) {
            addterm <- get_W_cdr_3D(arrayC(beta.new[[1]],dim=c(p[1],R[2])),arrayC(beta.new[[2]],dim=c(p[2],R[2])),arrayC(beta.new[[3]],dim=c(p[3],R[2])),
                                    tau.new, alpha.new)
          } else if (D==2) {
            addterm <- get_W_cdr_2D(arrayC(beta.new[[1]],dim=c(p[1],R[2])),arrayC(beta.new[[2]],dim=c(p[2],R[2])),
                                    tau.new, alpha.new)
          }
          addterm[which(addterm<.Machine$double.xmin)] <- .Machine$double.xmin
          c.dr_array_Bi <- c.dr_array_Bi + addterm
        }
        
        for (d in 1:D) {
          for (r in 1:R[2]) {
            #w.store[[2]][d,r,iter+1] <- max(rgig(1,1-.5*p[d]*Ni, c.dr_array_Bi[d,r], lambda.store[[2]][d,r,iter]),.Machine$double.xmin)
            w.store[[2]][d,r,iter+1] <- max(.Call("rgig",1,1-.5*p[d]*Ni, c.dr_array_Bi[d,r], lambda.store[[2]][d,r,iter],PACKAGE="GIGrvg"),.Machine$double.xmin)
            lambda.store[[2]][d,r,iter+1] <- rgamma(1,a.lam+p[d],b.lam+.5*p[d]*w.store[[2]][d,r,iter+1])
          }
        }

        alpha.cand <- arrayC(exp(rnorm(D*R[2], log(alpha.store[[2]][,,iter]), sigma_log_alpha)),dim=c(D,R[2]))
        acc_ratio <- exp(log_g_alpha_nvars2(alpha.cand, p, lapply(beta.store[[2]], function(x) x[,,,iter+1]),
                                            w.store[[2]][,,iter+1], tau.store[[2]][,iter+1], a.alpha, b.alpha, D, R[2], Ni)-
                           log_g_alpha_nvars2(arrayC(alpha.store[[2]][,,iter],dim=c(D,R[2])), p, lapply(beta.store[[2]], function(x) x[,,,iter+1]),
                                              w.store[[2]][,,iter+1], tau.store[[2]][,iter+1], a.alpha, b.alpha, D, R[2], Ni))
        alpha.store[[2]][,,iter+1] <- ifelse(arrayC(runif(D*R[2]),dim=c(D,R[2]))>acc_ratio, alpha.store[[2]][,,iter], alpha.cand)


      }
      
      # 3. update Gamma terms
      if (!null_Gamma) {
        if (show.prog && iter%%prog.count==0 && show.allsteps) {
          print(paste0(iter,': Updating Gamma'))
        }
        
        Y_G_vec <- getResid_fullModel(Yvec,ID,Time,Visit,Ci,Xi,Zti,
                                      M_iter,Bi_iter,0*Gamma_iter,Thetai_iter,Btm_iter_cond,
                                      Ds_iter,Gammak_iter,M)
        Y_G <- arrayC(Y_G_vec, dim=c(N,p))
        g(beta.new, w.new, lambda.new, tau.new, sigma2.new, alpha.new) %=% btrr.1var.1visit.1iter_rcpp(Y = Y_G, xc = as.matrix(Time), R=R[3], 
                                                                                                          beta.prev = lapply(beta.store[[3]],function(x) arrayC(x[,,iter],dim=c(dim(x)[1],R[3]))), w.prev=arrayC(w.store[[3]][,,iter],dim=c(D,R[3])), 
                                                                                                          lambda.prev = arrayC(lambda.store[[3]][,,iter],dim=c(D,R[3])), tau.prev = tau.store[[3]][iter], 
                                                                                                          sigma2.prev = sigma2.store[,iter], alpha.prev = arrayC(alpha.store[[3]][,,iter],dim=c(D,R[3])),
                                                                                                          update.w = T, update.lambda = T, update.alpha = T, update.tau = T, update.sigma2 = F, sigma_log_alpha = sigma_log_alpha)
        for (d in 1:D) {beta.store[[3]][[d]][,,iter+1] <- beta.new[[d]]}
        w.store[[3]][,,iter+1] <- w.new
        lambda.store[[3]][,,iter+1] <- lambda.new
        tau.store[[3]][iter+1] <- tau.new
        alpha.store[[3]][,,iter+1] <- alpha.new
        
        Gamma_iter <- c(TP.rankR(beta.new))
        
      }
      
      # 4. update Thetai terms
      if (!null_Thetai) {
        if (show.prog && iter%%prog.count==0 && show.allsteps) {
          print(paste0(iter,': Updating Thetai'))
        }
        c.dr_array_Thetai <- array(0,dim=c(D,R[4]))
        for (i in 1:Ni){
          Visit_i <- Visit[ID==i]
          Time_i <- Time[ID==i]
          Zti_i <- as.matrix(Zti[ID==i,])
          ID_i <- ID[ID==i]
          Y_Thetai_vec <- getResid_fullModel(Yvec[ID==i,],ID_i,Time_i,Visit_i,Ci,Xi,Zti_i,
                                             M_iter,Bi_iter,Gamma_iter,0*Thetai_iter,Btm_iter_cond,Ds_iter,Gammak_iter,M)
          Y_Thetai <- arrayC(Y_Thetai_vec, dim=c(dim(Y_Thetai_vec)[1],p))
          g(beta.new, w.new, lambda.new, tau.new, sigma2.new, alpha.new) %=% btrr.1var.1visit.1iter_rcpp(Y = Y_Thetai, xc = as.matrix(Time_i), R=R[4],
                                                                                                               beta.prev = lapply(beta.store[[4]],function(x) arrayC(x[,,i,iter],dim=c(dim(x)[1],R[4]))), w.prev=arrayC(w.store[[4]][,,iter],dim=c(D,R[4])),
                                                                                                               lambda.prev = arrayC(lambda.store[[4]][,,iter],dim=c(D,R[4])), tau.prev = tau.store[[4]][i,iter],
                                                                                                               sigma2.prev = sigma2.store[ID==i,iter], alpha.prev = arrayC(alpha.store[[4]][,,iter],dim=c(D,R[4])),
                                                                                                               update.w=F, update.lambda = F, update.tau = T, update.sigma2 = F, update.alpha = F)
          for (d in 1:D) {beta.store[[4]][[d]][,,i,iter+1] <- beta.new[[d]]}
          tau.store[[4]][i,iter+1] <- tau.new

          Thetai_iter[i,] <- c(TP.rankR(beta.new))
          
          if (D==3) {
            addterm <- get_W_cdr_3D(arrayC(beta.new[[1]],dim=c(p[1],R[4])),arrayC(beta.new[[2]],dim=c(p[2],R[4])),arrayC(beta.new[[3]],dim=c(p[3],R[4])),
                                    tau.new, alpha.new)
          } else if (D==2) {
            addterm <- get_W_cdr_2D(arrayC(beta.new[[1]],dim=c(p[1],R[4])),arrayC(beta.new[[2]],dim=c(p[2],R[4])),
                                    tau.new, alpha.new)
          }
          addterm[which(addterm<.Machine$double.xmin)] <- .Machine$double.xmin
          c.dr_array_Thetai <- c.dr_array_Thetai + addterm
        }
        
        # update w, lambda, and alpha
        for (d in 1:D) {
          for (r in 1:R[4]) {
            #w.store[[4]][d,r,iter+1] <- max(rgig(1, 1-.5*p[d]*Ni, c.dr_array_Thetai[d,r], lambda.store[[4]][d,r,iter]),.Machine$double.xmin)
            w.store[[4]][d,r,iter+1] <- max(.Call("rgig",1, 1-.5*p[d]*Ni, c.dr_array_Thetai[d,r], lambda.store[[4]][d,r,iter],PACKAGE="GIGrvg"),.Machine$double.xmin)
            lambda.store[[4]][d,r,iter+1] <- rgamma(1, a.lam+p[d], b.lam+.5*p[d]*w.store[[4]][d,r,iter+1])
          }
        }

        alpha.cand <- arrayC(exp(rnorm(D*R[4], log(alpha.store[[4]][,,iter]), sigma_log_alpha)),dim=c(D,R[4]))
        acc_ratio <- exp(log_g_alpha_nvars2(alpha.cand, p, lapply(beta.store[[4]], function(x) x[,,,iter+1]),
                                            w.store[[4]][,,iter+1], tau.store[[4]][,iter+1], a.alpha, b.alpha, D, R[4], Ni)-
                           log_g_alpha_nvars2(arrayC(alpha.store[[4]][,,iter],dim=c(D,R[4])), p, lapply(beta.store[[4]], function(x) x[,,,iter+1]),
                                              w.store[[4]][,,iter+1], tau.store[[4]][,iter+1], a.alpha, b.alpha, D, R[4], Ni))
        alpha.store[[4]][,,iter+1] <- ifelse(arrayC(runif(D*R[4]),dim=c(D,R[4]))>acc_ratio, alpha.store[[4]][,,iter], alpha.cand)

      }
      
      
      # 5. update B_{tm} terms
      if (!null_Btm) {
        if (show.prog && iter%%prog.count==0 && show.allsteps) {
          print(paste0(iter,': Updating B_tm'))
        }
        
        c.dr_array_Btm <- array(0,dim=c(D,R[5]))
        for (t in ifelse(null_baseline,2,1):Nt) {
          ID_t <- ID[Visit==(t-1)]
          Zti_t <- as.matrix(Zti[Visit==(t-1),])
          Time_t <- Time[Visit==(t-1)]
          Visit_t <- Visit[Visit==(t-1)]
          for (m in 1:M) {
            Btm_iter_mstar <- Btm_iter
            Btm_iter_mstar[t,m,] <- 0
            Y_Btm_vec <- getResid_fullModel(Yvec[Visit==(t-1),],ID_t,Time_t,Visit_t,Ci,Xi,Zti_t,
                                            M_iter,Bi_iter,Gamma_iter,Thetai_iter,arrayC(Btm_iter_mstar,dim=c(Nt*M,prod(p))),
                                            Ds_iter,Gammak_iter,M)
            Y_Btm <- arrayC(Y_Btm_vec, dim=c(dim(Y_Btm_vec)[1],p))
            g(beta.new, w.new, lambda.new, tau.new, sigma2.new, alpha.new) %=% btrr.1var.1visit.1iter_rcpp(Y=Y_Btm, xc=as.matrix(Ci[ID_t,m]), R=R[5],
                                                                                                              beta.prev=lapply(beta.store[[5]],function(x) arrayC(x[,,t,m,iter],dim=c(dim(x)[1],R[5]))), w.prev=arrayC(w.store[[5]][,,iter],dim=c(D,R[5])),
                                                                                                              lambda.prev=arrayC(lambda.store[[5]][,,iter],dim=c(D,R[5])),tau.prev=tau.store[[5]][t,m,iter],
                                                                                                              sigma2.prev=sigma2.store[Visit==(t-1),iter], alpha.prev=arrayC(alpha.store[[5]][,,iter],dim=c(D,R[5])),
                                                                                                              update.w = F, update.lambda = F, update.tau = T, update.sigma2 = F, update.alpha = F)
            for (d in 1:D) {beta.store[[5]][[d]][,,t,m,iter+1] <- beta.new[[d]]}
            tau.store[[5]][t,m,iter+1] <- tau.new
            
            Btm_iter[t,m,] <- c(TP.rankR(beta.new))
            
            if (D==3) {
              addterm <- get_W_cdr_3D(arrayC(beta.new[[1]],dim=c(p[1],R[5])),arrayC(beta.new[[2]],dim=c(p[2],R[5])),arrayC(beta.new[[3]],dim=c(p[3],R[5])),
                                      tau.new, alpha.new)
            } else if (D==2) {
              addterm <- get_W_cdr_2D(arrayC(beta.new[[1]],dim=c(p[1],R[5])),arrayC(beta.new[[2]],dim=c(p[2],R[5])),
                                      tau.new, alpha.new)
            }
            addterm[which(addterm<.Machine$double.xmin)] <- .Machine$double.xmin
            c.dr_array_Btm <- c.dr_array_Btm+addterm
            
          }
        }
        Btm_iter_cond <- arrayC(Btm_iter,dim=c(Nt*M,prod(p)))
        
        # update w, lambda, and alpha
        for (d in 1:D) {
          for (r in 1:R[5]) {
            #w.store[[5]][d,r,iter+1] <- max(rgig(1, 1-.5*p[d]*Nt*M, c.dr_array_Btm[d,r], lambda.store[[5]][d,r,iter]),.Machine$double.xmin)
            w.store[[5]][d,r,iter+1] <- max(.Call("rgig",1, 1-.5*p[d]*Nt*M, c.dr_array_Btm[d,r], lambda.store[[5]][d,r,iter],PACKAGE="GIGrvg"),.Machine$double.xmin)
            lambda.store[[5]][d,r,iter+1] <- rgamma(1, a.lam+p[d], b.lam+.5*p[d]*w.store[[5]][d,r,iter+1])
          }
        }
        
        alpha.cand <- arrayC(exp(rnorm(D*R[5], log(alpha.store[[5]][,,iter]), sigma_log_alpha)),dim=c(D,R[5]))
        acc_ratio <- exp(log_g_alpha_nvars_nvisits2(alpha.cand, p, lapply(beta.store[[5]], function(x) x[,,,,iter+1]),
                                                    w.store[[5]][,,iter+1], tau.store[[5]][,,iter+1], a.alpha, b.alpha, D, R[5], Nt, M)-
                           log_g_alpha_nvars_nvisits2(arrayC(alpha.store[[5]][,,iter],dim=c(D,R[5])), p, lapply(beta.store[[5]], function(x) x[,,,,iter+1]),
                                                      w.store[[5]][,,iter+1], tau.store[[5]][,,iter+1], a.alpha, b.alpha, D, R[5], Nt,M))
        alpha.store[[5]][,,iter+1] <- ifelse(arrayC(runif(D*R[5]),dim=c(D,R[5]))>acc_ratio, alpha.store[[5]][,,iter], alpha.cand)
        
      }
      
      # 6. update D_s terms
      if (!null_Ds) {
        if (show.prog && iter%%prog.count==0 && show.allsteps) {
          print(paste0(iter,': Updating D_s'))
        }
        c.dr_array_Ds <- array(0,dim=c(D,R[6]))
        for (s in 1:S) {
          Ds_iter_sstar <- Ds_iter
          Ds_iter_sstar[s,] <- 0
          Y_Ds_vec <- getResid_fullModel(Yvec,ID,Time,Visit,Ci,Xi,Zti,
                                         M_iter,Bi_iter,Gamma_iter,Thetai_iter,Btm_iter_cond,
                                         Ds_iter_sstar,Gammak_iter,M)
          Y_Ds <- arrayC(Y_Ds_vec, dim=c(N,p))
          g(beta.new, w.new, lambda.new, tau.new, sigma2.new, alpha.new) %=% btrr.1var.1visit.1iter_rcpp(Y=Y_Ds, xc=as.matrix(sapply(ID, function(x) Xi[x,s])), R=R[6],
                                                                                                            beta.prev=lapply(beta.store[[6]],function(x) arrayC(x[,,s,iter],dim=c(dim(x)[1],R[6]))), w.prev=arrayC(w.store[[6]][,,iter],dim=c(D,R[6])),
                                                                                                            lambda.prev=arrayC(lambda.store[[6]][,,iter],dim=c(D,R[6])),tau.prev=tau.store[[6]][s,iter],
                                                                                                            sigma2.prev=sigma2.store[,iter], alpha.prev=arrayC(alpha.store[[6]][,,iter],dim=c(D,R[6])),
                                                                                                            update.w = F, update.lambda = F, update.tau = T, update.sigma2 = F, update.alpha = F)
          for (d in 1:D) {beta.store[[6]][[d]][,,s,iter+1] <- beta.new[[d]]}
          tau.store[[6]][s,iter+1] <- tau.new
          
          Ds_iter[s,] <- c(TP.rankR(beta.new))
          
          if (D==3) {
            addterm <- get_W_cdr_3D(arrayC(beta.new[[1]],dim=c(p[1],R[6])),arrayC(beta.new[[2]],dim=c(p[2],R[6])),arrayC(beta.new[[3]],dim=c(p[3],R[6])),
                                    tau.new, alpha.new)
          } else if (D==2) {
            addterm <- get_W_cdr_2D(arrayC(beta.new[[1]],dim=c(p[1],R[6])),arrayC(beta.new[[2]],dim=c(p[2],R[6])),
                                    tau.new, alpha.new)
          }
          addterm[which(addterm<.Machine$double.xmin)] <- .Machine$double.xmin
          c.dr_array_Ds <- c.dr_array_Ds + addterm
        }
        
        # update w, lambda, and alpha
        for (d in 1:D) {
          for (r in 1:R[6]) {
            #w.store[[6]][d,r,iter+1] <- max(rgig(1, 1-.5*p[d]*S, c.dr_array_Ds[d,r], lambda.store[[6]][d,r,iter]),.Machine$double.xmin)
            w.store[[6]][d,r,iter+1] <- max(.Call("rgig",1, 1-.5*p[d]*S, c.dr_array_Ds[d,r], lambda.store[[6]][d,r,iter],PACKAGE="GIGrvg"),.Machine$double.xmin)
            lambda.store[[6]][d,r,iter+1] <- rgamma(1, a.lam+p[d], b.lam+.5*p[d]*w.store[[6]][d,r,iter+1])
          }
        }
        
        alpha.cand <- arrayC(exp(rnorm(D*R[6], log(alpha.store[[6]][,,iter]), sigma_log_alpha)),dim=c(D,R[6]))
        acc_ratio <- exp(log_g_alpha_nvars2(alpha = alpha.cand, p = p, beta = lapply(beta.store[[6]], function(x) x[,,,iter+1]),
                                            w = w.store[[6]][,,iter+1], tau = tau.store[[6]][,iter+1], a.alpha = a.alpha, b.alpha = b.alpha, D = D, R = R[6], K = S)-
                           log_g_alpha_nvars2(arrayC(alpha.store[[6]][,,iter],dim=c(D,R[6])), p, lapply(beta.store[[6]], function(x) x[,,,iter+1]),
                                              w.store[[6]][,,iter+1], tau.store[[6]][,iter+1], a.alpha, b.alpha, D, R[6], S))
        alpha.store[[6]][,,iter+1] <- ifelse(arrayC(runif(D*R[6]),dim=c(D,R[6]))>acc_ratio, alpha.store[[6]][,,iter], alpha.cand)
        
      }
      
      # 7. update Gamma_k terms
      if (!null_Gammak) {
        if (show.prog && iter%%prog.count==0 && show.allsteps) {
          print(paste0(iter,': Updating Gamma_k'))
        }
        c.dr_array_Gammak <- array(0,dim=c(D,R[7]))
        
        for (k in 1:K){
          Gammak_iter_kstar <- Gammak_iter
          Gammak_iter_kstar[k,] <- 0

          Y_Gk_vec <- getResid_fullModel(Yvec,ID,Time,Visit,Ci,Xi,Zti,
                                             M_iter,Bi_iter,Gamma_iter,Thetai_iter,Btm_iter_cond,Ds_iter,Gammak_iter_kstar,M)
          Y_Gk <- arrayC(Y_Gk_vec, dim=c(N,p))
          g(beta.new, w.new, lambda.new, tau.new, sigma2.new, alpha.new) %=% btrr.1var.1visit.1iter_rcpp(Y = Y_Gk, xc = as.matrix(Zti[,k]), R=R[7], 
                                                                                                               beta.prev = lapply(beta.store[[7]],function(x) arrayC(x[,,k,iter],dim=c(dim(x)[1],R[7]))), w.prev=arrayC(w.store[[7]][,,iter],dim=c(D,R[7])),
                                                                                                               lambda.prev = arrayC(lambda.store[[7]][,,iter],dim=c(D,R[7])), tau.prev = tau.store[[7]][k,iter],
                                                                                                               sigma2.prev = sigma2.store[,iter], alpha.prev = arrayC(alpha.store[[7]][,,iter],dim=c(D,R[7])),
                                                                                                               update.w = F, update.lambda = F, update.tau = T, update.sigma2 = F, update.alpha = F)
          for (d in 1:D) {beta.store[[7]][[d]][,,k,iter+1] <- beta.new[[d]]}
          tau.store[[7]][k,iter+1] <- tau.new
          
          Gammak_iter[k,] <- c(TP.rankR(beta.new))
          
          if (D==3) {
            addterm <- get_W_cdr_3D(arrayC(beta.new[[1]],dim=c(p[1],R[7])),arrayC(beta.new[[2]],dim=c(p[2],R[7])),arrayC(beta.new[[3]],dim=c(p[3],R[7])),
                                    tau.new, alpha.new)
          } else if (D==2) {
            addterm <- get_W_cdr_2D(arrayC(beta.new[[1]],dim=c(p[1],R[7])),arrayC(beta.new[[2]],dim=c(p[2],R[7])),
                                    tau.new, alpha.new)
          }
          addterm[which(addterm<.Machine$double.xmin)] <- .Machine$double.xmin
          c.dr_array_Gammak <- c.dr_array_Gammak + addterm
        }
        # update w, lambda, and alpha
        for (d in 1:D) {
          for (r in 1:R[7]) {
            #w.store[[7]][d,r,iter+1] <- max(rgig(1, 1-.5*p[d]*K, c.dr_array_Gammak[d,r], lambda.store[[7]][d,r,iter]),.Machine$double.xmin)
            w.store[[7]][d,r,iter+1] <- max(.Call("rgig",1, 1-.5*p[d]*K, c.dr_array_Gammak[d,r], lambda.store[[7]][d,r,iter],PACKAGE="GIGrvg"),.Machine$double.xmin)
            lambda.store[[7]][d,r,iter+1] <- rgamma(1, a.lam+p[d], b.lam+.5*p[d]*w.store[[7]][d,r,iter+1])
          }
        }
        
        alpha.cand <- arrayC(exp(rnorm(D*R[7], log(alpha.store[[7]][,,iter]), sigma_log_alpha)),dim=c(D,R[7]))
        acc_ratio <- exp(log_g_alpha_nvars2(alpha.cand, p, lapply(beta.store[[7]], function(x) x[,,,iter+1]),
                                            w.store[[7]][,,iter+1], tau.store[[7]][,iter+1], a.alpha, b.alpha, D, R[7], K)-
                           log_g_alpha_nvars2(arrayC(alpha.store[[7]][,,iter],dim=c(D,R[7])), p, lapply(beta.store[[7]], function(x) x[,,,iter+1]),
                                              w.store[[7]][,,iter+1], tau.store[[7]][,iter+1], a.alpha, b.alpha, D, R[7], K))
        alpha.store[[7]][,,iter+1] <- ifelse(arrayC(runif(D*R[7]),dim=c(D,R[7]))>acc_ratio, alpha.store[[7]][,,iter], alpha.cand)
        
      }

      # update sigma^2
      if (show.prog && iter%%prog.count==0 && show.allsteps) {
        print(paste0(iter, ': Updating sigma^2'))
      }
      sigma2.store[,iter+1] <- sampleSigma_fullModel(Yvec,ID,Time,Visit,Ci,Xi,Zti,M_iter,Bi_iter,
                                                     Gamma_iter,Thetai_iter,Btm_iter_cond,Ds_iter,Gammak_iter,M,a.sig,b.sig)
    }
  return(list(beta.store, w.store, lambda.store, tau.store, sigma2.store, alpha.store))
  
}
