#library(tidyverse)
library(dplyr)
library(tidyr)
library(gsubfn)
#library(plot.matrix)
library(MASS)
library(tictoc)
library(GIGrvg)
library(rTensor)
library(lmerTest)
library(rlist)

### MAIN FUNCTIONS START HERE

# tensor.reg: Full tensor regression over nsweep iterations using rank-R decomposition; includes random intercept
# Estimates models of the form: Y_i = \Gamma_1 x1 + ... + \Gamma_K xK + Bi + \epsilon, where Y_i is a tensor outcome, 
# x1,...,xK are scalar predictors, \Gamma's are tensor coefficients, and Bi is a subject-specific intercept
#
# INPUTS:
#     Y.train -- Array containing tensor outcome information. 
#                Dimension: N x p1 x ... x pD, where N is sample size and p1-pD are tensor dimension sizes up to dimension D
#     X.train -- Array containing scalar covariate information.
#                Dimension: N x K, where N is sample size and K is number of covariates
#     ID      -- Vector of integers of length N. Values must range from 1 to number of unique subjects. 
#                Set ID=1 or ID=NA to not estimate subject-specific intercepts
#     nsweep  -- Number of MCMC iterations to be performed
#     R       -- Rank of tensor decomposition placed on tensor coefficients. Higher ranks have better accuracy, but longer run-times.
#                If subject-specific intercepts are estimated, you may specify rank of those terms with second element of R, e.g. R=c(3,2).
#                Otherwise, rank of main terms and subject-specific intercepts are assumed equal.
#     a.sig/b.sig -- Hyper-parameters for noise term (Inverse Gamma prior)
#     a.tau/b.tau -- Hyper-parameters for tensor margin variance scaling term (Gamma prior)
#     a.lam/b.lam -- Hyper-parameters for rate parameter (Gamma prior)
#     init.type   -- Set to "ols" if OLS should be used for initialization of model parameters. By default, does not perform OLS.
#                    Currently only works when no subject-specific intercepts are included (i.e. ID=NA)
#     show.prog   -- Set to T if number of iterations should be periodically displayed as MCMC is performed. Default is T.
#     prog.count  -- Number specifying how often to display MCMC iteration count if show.prog=T. Default is 10.
#
# OUTPUTS: results <- tensor.reg(Y.train, X.train, ...)
#     Overall Parameters:
#     results[[1]] -- Array containing all tensor regression coefficients for each variable, at each iteration of MCMC
#                     Dimension: nsweep x K x (p1...pD). 
#                     e.g. results[[1]][5,1,] is vectorized version of 1st tensor coefficient at the 5th iteration of MCMC
#                     To obtain the mean estimate over all iterations, use apply(results[[1]], c(2,3), mean).
#                     To obtain credible intervals, use apply(results[[1]], c(2,3), function(x) quantile(c(lower,upper)))
#     results[[2]] -- List of length D (number of dimensions), each element contains an array of tensor margins
#                     dim(results[[2]][[d]]) is nsweep x pd x R x K, where pd is size of d-th dimension and K is number of variables
#     results[[3]] -- List of length D, each element contains an array of tensor margin covariances
#                     dim(results[[3]][[d]]) is nsweep x pd x R x K
#     results[[4]] -- Array containing lambda terms (rate parameter). Dimension: nsweep x D x R x K, where D is number of dimensions
#     results[[5]] -- Array containing tau terms (variance scaling). Dimension: nsweep x K, where K is number of variables
#     results[[6]] -- Vector containing sigma terms (noise variance) across all iterations. Length: nsweep.
#                     Useful for assessing convergence in traceplots, i.e. plot(results[[6]]).

#     Subject-specific Parameters (NA if ID=NA/not specified)
#     results[[7]] -- Array containing all vectorized ID-specific intercepts (tensors), at each iteration of MCMC
#                     Dimension: nsweep x n_IDs x (p1...pD), where n_IDs is number of unique subjects
#     results[[8]] -- List of length D, each element contains an array of tensor margins for the subject-specific intercepts
#                     dim(results[8]][[d]]) is nsweep x pd x R x n_IDs, where pd is size of d-th dimension and n_IDs is number of subjects
#     results[[9]] -- List of length D, each element contains an array of tensor margin covariances for the subject-specific intercepts
#                     dim(results[9]][[d]]) is nsweep x pd x R x n_IDs, where pd is size of d-th dimension and n_IDs is number of subjects
#     results[[10]] -- Array containing lambda terms (rate parameter) corresponding to subject-specific intercepts. 
#                     Dimension: nsweep x D x R x n_IDs, where n_IDs is number of subjects

tensor.reg <- function(Y.train, X.train, ID=NA, nsweep=1000, R=3,
                       a.sig=1, b.sig=1, a.tau=1, b.tau=1,
                       a.lam, b.lam, init.type, show.prog=T, prog.count=10) {
  # require necessary packages
  require(GIGrvg)
  
  # set size parameters and defaults
  N <- dim(Y.train)[1]
  p <- dim(Y.train)[-1]
  D <- length(p)
  K <- ifelse(is.null(ncol(X.train)), 1, ncol(X.train))
  ngroups <- length(unique(ID))
  if (missing(a.lam)) {a.lam <- D}
  if (missing(b.lam)) {b.lam <- a.lam ** (1/(2*D))}
  if (missing(init.type)) {init.type <- "default"}
  if (length(R)==1) {R2 <- R} # default is same
  if (length(R)>1) {R2 <- R[2]; R <- R[1]} # rank specification for subject-specific tensor coefs
  
  # initialize model parameters
  margins.init <- vector(mode="list",D)
  margins.cov.init <- vector(mode="list",D)
  margins.ri.init <- vector(mode="list",D)
  margins.ri.cov.init <- vector(mode="list",D)
  for (d in 1:D) {
    var.d <- rgamma(1,a.tau,b.tau)
    margins.init[[d]] <- array(rnorm(p[d]*R*K,sd=sqrt(var.d)), dim=c(p[d],R,K))
    margins.cov.init[[d]] <- array(1,dim=c(p[d],R,K))
    margins.ri.init[[d]] <- array(rnorm(p[d]*R2*ngroups),dim=c(p[d],R2,ngroups))
    margins.ri.cov.init[[d]] <- array(1,dim=c(p[d],R2,ngroups))
  }
  lambda.init <- array(1,dim=c(D,R,K))
  lambda.ri.init <- array(1,dim=c(D,R2,ngroups))
  tau.init <- rep(1,K)
  sigma2.init <- 1
  
  # use OLS to get initial margins if specified
  if (init.type == "ols") {
    # Run voxel-wise ols on Y.train ~ X.train , ignoring random effects
    Gamma.allk.ols <- vox.ols.reg(Y.train,X.train)
    
    # take tensor decomposition of each Gamma.k to get initial margins
    for (k in 1:K) {
      margins.k.ols <- tensor.decomp(array(Gamma.allk.ols[k,], dim=p),rank=R)
      for (d in 1:D){
        margins.init[[d]][,,k] <- margins.k.ols[[d]]
      }
    }
  }
  
  # get initial coefficient estimates (Gamma's)
  Gamma.all.vec.init <- array(NA,dim=c(K,prod(p)))
  for (k in 1:K) {
    Gamma.all.vec.init[k,] <- c(TP.rankR(lapply(margins.init, function(x) x[,,k])))
  }
  Gamma.ri.vec.init <- array(NA,dim=c(ngroups,prod(p)))
  for (i in 1:ngroups){
    Gamma.ri.vec.init[i,] <- c(TP.rankR(lapply(margins.ri.init, function(x) x[,,i])))
  }
  
  # initialize storage variables
  Gamma.allk.store <- array(NA,dim=c(nsweep,K,prod(p)))
  Gamma.ri.store <- array(NA,dim=c(nsweep,ngroups,prod(p))) # ri=random intercept
  margins.store <- vector(mode="list",D)
  margins.cov.store <- vector(mode="list",D)
  margins.ri.store <- vector(mode="list",D) #
  margins.ri.cov.store <- vector(mode="list",D) #
  for (d in 1:D) {
    margins.store[[d]] <- array(NA,dim=c(nsweep,p[d],R,K))
    margins.cov.store[[d]] <- array(NA,dim=c(nsweep,p[d],R,K))
    margins.ri.store[[d]] <- array(NA,dim=c(nsweep,p[d],R2,ngroups)) #
    margins.ri.cov.store[[d]] <- array(NA,dim=c(nsweep,p[d],R2,ngroups)) #
  }
  lambda.store <- array(NA,dim=c(nsweep,D,R,K))
  lambda.ri.store <- array(NA,dim=c(nsweep,D,R2,ngroups)) #
  tau.store <- array(NA,dim=c(nsweep,K))
  sigma2.store <- rep(NA,nsweep)
  
  # put initial parameters as first iteration of storage variables
  Gamma.allk.store[1,,] <- Gamma.all.vec.init
  Gamma.ri.store[1,,] <- Gamma.ri.vec.init
  for (d in 1:D) {
    margins.store[[d]][1,,,] <- margins.init[[d]]
    margins.cov.store[[d]][1,,,] <- margins.cov.init[[d]]
    margins.ri.store[[d]][1,,,] <- margins.ri.init[[d]]
    margins.ri.cov.store[[d]][1,,,] <- margins.ri.cov.init[[d]]
  }
  lambda.store[1,,,] <- lambda.init
  lambda.ri.store[1,,,] <- lambda.ri.init
  tau.store[1,] <- tau.init
  sigma2.store[1] <- sigma2.init
  
  
  # loop over iterations
  # note: treat initialization as iteration=1
  for (s in 2:nsweep) {
    
    if (show.prog && s%%prog.count==0) {print(s)}
    
    for (k in 1:K) {
      Y.tilde.k.vec <- array(Y.train,dim=c(N,prod(p)))
      for (k.star in 1:K) {
        # use previous iteration if k.star > k, and current if k.star < k
        if (k.star > k) {Gamma.k <- Gamma.allk.store[s-1,k.star,]}
        if (k.star < k) {Gamma.k <- Gamma.allk.store[s,k.star,]}
        # use previous iteration for all Gamma.k's
        # if (k != k.star) {Gamma.k <- Gamma.allk.store[s-1,k.star,]}
        if (k == k.star) {Gamma.k <- 0}
        for (n in 1:N) {
          Y.tilde.k.vec[n,] <- Y.tilde.k.vec[n,] - X.train[n,k.star]*c(Gamma.k)
        }
      }
      # remove subject-specific intercept if applicable
      if (ngroups > 1 && ngroups < N) {
        for (n in 1:N) {
          I <- ID[n]
          Y.tilde.k.vec[n,] <- Y.tilde.k.vec[n,] - Gamma.ri.store[s-1,I,]
        }
      }
      
      # run STR to get Gamma_k
      Y.tilde.k <- array(Y.tilde.k.vec, dim=c(N,p))
      xk.train <- X.train[,k]
      margins.k <- lapply(margins.store, function(x) x[s-1,,,k])
      margins.cov.k <- lapply(margins.cov.store, function(x) x[s-1,,,k])
      lambda.k <- lambda.store[s-1,,,k]
      tau.k <- tau.store[s-1,k]
      sigma2 <- sigma2.store[s-1]
      
      STR_1iter_results <- tensor.reg.1var.1iter(Y.tilde.k, xk.train, R=R, a.sig=a.sig, b.sig=b.sig,
                                                 a.tau=a.tau, b.tau=b.tau, a.lam=a.lam, b.lam=b.lam,
                                                 margins.prev=margins.k, margins.cov.prev=margins.cov.k,
                                                 lambda.prev=lambda.k, tau.prev=tau.k,
                                                 sigma2.prev=sigma2)
      
      # store k-specific parameters
      Gamma.allk.store[s,k,] <- c(STR_1iter_results[[1]])
      for (d in 1:D) {
        margins.store[[d]][s,,,k] <- STR_1iter_results[[2]][[d]]
        margins.cov.store[[d]][s,,,k] <- STR_1iter_results[[3]][[d]]
      }
      lambda.store[s,,,k] <- STR_1iter_results[[4]]
      tau.store[s,k] <- STR_1iter_results[[5]]
    }
    
    # if 1 < ngroups < N, sample subject-specific intercept terms by dividing data into groups
    if (ngroups > 1 && ngroups < N) {
      for (i in 1:ngroups) {
        Ni <- sum(ID==i)
        Yi.vec <- array(Y.train, dim=c(N,prod(p)))[ID==i,]
        # remove all main effects
        for (k in 1:K) {
          Gamma.k <- Gamma.allk.store[s,k,]
          for (n in Ni) {
            Yi.vec[n,] <- Yi.vec[n,] - X.train[n,k]*Gamma.k
          }
        }
        # run STR to find intercept for i-th group after removing other effects
        Y.ri <- array(Yi.vec, dim=c(Ni,p))
        x.ri <- rep(1,Ni)
        margins.ri <- lapply(margins.ri.store, function(x) x[s-1,,,i])
        margins.cov.ri <- lapply(margins.ri.cov.store, function(x) x[s-1,,,i])
        lambda.ri <- lambda.ri.store[s-1,,,i]
        tau.ri <- 1
        sigma2.ri <- sigma2.store[s-1]
        
        STR_ri_results <- tensor.reg.1var.1iter(Y.ri, x.ri, R=R2, a.sig=a.sig, b.sig=b.sig,
                                                a.tau=a.tau, b.tau=b.tau, a.lam=a.lam, b.lam=b.lam,
                                                margins.prev=margins.ri, margins.cov.prev=margins.cov.ri,
                                                lambda.prev=lambda.ri,tau.prev=tau.ri,sigma2.prev=sigma2.ri)
        # store random intercept parameters
        Gamma.ri.store[s,i,] <- c(STR_ri_results[[1]])
        for (d in 1:D) {
          margins.ri.store[[d]][s,,,i] <- STR_ri_results[[2]][[d]]
          margins.ri.cov.store[[d]][s,,,i] <- STR_ri_results[[3]][[d]]
        }
        lambda.ri.store[s,,,i] <- STR_ri_results[[4]]
      }
    }
    
    # find residual to update sigma2
    Y.train.vec <- array(Y.train, dim=c(N,prod(p)))
    Resid.vec <- Y.train.vec
    tot_vox <- prod(dim(Y.train.vec)) - sum(is.na(Y.train.vec))
    
  
    for (n in 1:N) {
      # remove main effects
      for (k in 1:K) {
        Resid.vec[n,] <- Resid.vec[n,] - X.train[n,k]*Gamma.allk.store[s,k,] 
      }
      # remove id-specific intercepts
      if (ngroups > 1 && ngroups < N) {
        id <- ID[n]
        Resid.vec[n,] <- Resid.vec[n,] - Gamma.ri.store[s,id,]
      }
    }
    
    # missing values have no contribution towards posterior
    sigma2.store[s] <- 1/rgamma(1, a.sig + tot_vox / 2, 
                                b.sig + .5 * sum(Resid.vec^2, na.rm=T))
    
  }
  
  return(list(Gamma.allk.store, margins.store, margins.cov.store, lambda.store,
              tau.store, sigma2.store, Gamma.ri.store, margins.ri.store,
              margins.ri.cov.store, lambda.ri.store))
  
}


# tensor.reg.1var.1iter: Single iteration for 1-variable tensor regression case which can accommodate missing outcome values
tensor.reg.1var.1iter <- function(Y.train, x.train, R=3,
                                  a.sig=1, b.sig=1, a.tau=1, b.tau=1,
                                  a.lam, b.lam, margins.prev, margins.cov.prev,
                                  lambda.prev, tau.prev, sigma2.prev) {
  
  N <- length(x.train)
  p <- dim(Y.train)[-1]
  D <- length(p)
  
  # set default hyper-parameters
  if (missing(a.lam)) {a.lam <- D}
  if (missing(b.lam)) {b.lam <- a.lam ** (1/(2*D))}
  
  # If R=1, set margins to be matrices if not already
  if (R==1) {
    margins.prev <- lapply(margins.prev, function(x) as.matrix(x))
    margins.cov.prev <- lapply(margins.cov.prev, function(x) as.matrix(x))
    lambda.prev <- as.matrix(lambda.prev)
  }
  
  # initialized updated values
  margins.new <- margins.prev
  margins.cov.new <- margins.cov.prev
  lambda.new <- lambda.prev
  tau.new <- tau.prev
  sigma2.new <- sigma2.prev
  
  norm.term.tau <- 0
  
  # update rank-specific terms
  for (r in 1:R) {
    
    # account for changes in previous r terms by redefining Gamma.allr
    Gamma.allr <- TP.rankR(margins.new)
    Gamma.rother <- 0*Gamma.allr
    if (R>1) {
      Gamma.rother <- Gamma.allr - TP.list(lapply(margins.new, function(x) x[,r]))
    }
    
    # remove other rank effects from Y.train
    Y.train.vec <- array(Y.train,dim=c(N,prod(p)))
    Y.tilde.vec <- array(NA,dim=c(N,prod(p)))
    for (n in 1:N) {
      Y.tilde.vec[n,] <- Y.train.vec[n,] - x.train[n] * c(Gamma.rother)
    }
    Y.tilde <- array(Y.tilde.vec, dim=c(N,p))
    
    
    # For each margin (d), update tensor margin (gamma), margin covariance (W), and
    # rate parameter (lambda)
    
    for (d in 1:D) {
      A_dr <- rep(0,p[d]) # defines posterior mean and covariance
      B_dr <- rep(0,p[d]) # defines posterior mean
      C_dr <- TP.list(lapply(margins.new[-d], function(x) x[,r]))
      
      for (i in 1:p[d]) {
        Y_di <- index.tensor(Y.tilde,d+1,i) # fix dimension d at element i
        A_dr[i] <- sigma2.new/(tau.new*margins.cov.prev[[d]][i,r])
        for (n in 1:N) {
          # convert missing indices to 0s (no contribution to C_drni)
          C_drni <- C_dr
          C_drni[is.na(index.tensor(Y_di,1,n))] <- 0
          Y_di_NAto0 <- Y_di %>% replace(is.na(.),0)
          #Y_di_NAto0 <- Y_di; Y_di_NAto0[is.na(Y_di_NAto0)] <- 0
          # add the term xn^2*sum(Cdni^2) to Ad[i]
          A_dr[i] <- A_dr[i] + x.train[n]^2*sum(C_drni^2)
          # add the term xn*(sum_{\Omega_v}(Yn*Cdni))
          B_dr[i] <- B_dr[i] + x.train[n]*sum(index.tensor(Y_di_NAto0,1,n)*C_drni)
        }
      }
      
      mu_dr <- B_dr / A_dr
      Sig_dr <- sigma2.new / A_dr
      
      # gamma_dr step: Draw from Normal
      gamma_dr <- rnorm(p[d], mean=mu_dr, sd=sqrt(Sig_dr))
      margins.new[[d]][,r] <- gamma_dr
      
      # W_dr step: Draw from generalized inverse Gaussian
      W_dr <- rgig(n=p[d], 1/2, gamma_dr^2 / tau.new, lambda.new[d,r])
      margins.cov.new[[d]][,r] <- W_dr
      
      # lambda_dr step: Draw from Gamma
      lambda_dr <- rgamma(1, a.lam + p[d], rate=b.lam+.5*sum(W_dr))
      lambda.new[d,r] <- lambda_dr
      
      # Add to norm term for updating tau
      norm.term.tau <- norm.term.tau + sum(gamma_dr^2 / W_dr)
    }
  }
  
  # tau step: Draw from generalized inverse Gaussian
  tau.new <- rgig(1, a.tau-.5*R*sum(p), norm.term.tau, 2*b.tau)
  
  
  # Update sigma
  Gamma.new <- TP.rankR(margins.new)
  tot_vox <- 0
  norm.term.sig <- 0
  for (n in c(1:N)) {
    tot_vox <- tot_vox + sum(!is.na(index.tensor(Y.train,1,n)))
    Yn.NAto0 <- index.tensor(Y.train,1,n) %>% replace(is.na(.),0)
    Gamma.n.NAto0 <- Gamma.new
    Gamma.n.NAto0[Yn.NAto0==0] <- 0
    norm.term.sig <- norm.term.sig + sum( (Yn.NAto0-x.train[n]*Gamma.n.NAto0)^2)
  }
  
  # sigma2 step: Draw from inverse Gamma
  sigma2.new <- 1/rgamma(1, a.sig+(tot_vox)/2, rate=b.sig+norm.term.sig/2)
  
  return(list(Gamma.new, margins.new, margins.cov.new, lambda.new, tau.new, sigma2.new))
}

### HELPER FUNCTIONS START HERE

# Given a list of arrays containing each dimension's tensor margins across all ranks (dim(X.allr[[d]]) = pd x R),
# returns the low-rank tensor made by summing up tensor products for all ranks
TP.rankR <- function(X.allr) { 
  R <- ncol(X.allr[[1]])
  if (is.null(R)) {
    return(TP.list(X.allr))
  } else {
    Y <- array(0, dim=c(as.numeric(lapply(X.allr, function(x) length(x[,1])))))
    for (r in c(1:R)) {
      Y <- Y + TP.list(lapply(X.allr, function(x) x[,r]))
    }
    return(Y)
  }
}

# Given a list of tensor margins (vectors), returns tensor product
TP.list <- function(X) {
  D <- length(X)
  if (D==1) {
    return(X[[1]])
  }
  if (D==2) {
    return(outer(X[[1]], X[[2]]))
  }
  else {
    return(outer(TP.list(X[1:(D-1)]), X[[D]]))
  }
}

# Fix a particular margin of a tensor X at a certain index
index.tensor <- function(X, margin, index) {
  p <- dim(X)
  D <- length(p)
  # code smaller cases so that run-time isn't as long
  if (D==1) {
    return(X[index])
  } else if (D==2) {
    if (margin==1) {
      return(X[index,])
    } else {
      return(X[,index])
    }
  } else if (D==3) {
    if (margin==1) {
      return(X[index,,])
    } else if (margin==2) {
      return(X[,index,])
    } else {
      return(X[,,index])
    }
  } else if (D==4) {
    if (margin==1) {
      return(X[index,,,])
    } else if (margin==2) {
      return(X[,index,,])
    } else if (margin==3) {
      return(X[,,index,])
    } else {
      return(X[,,,index])
    }
  # D > 4 case
  } else {
    X_dimnames <- vector(mode="list", length=D)
    for (d in c(1:D)) {
      X_dimnames[[d]] <- paste0(d,1:p[d])
    }
    X <- array(X, dim=dim(X), dimnames=X_dimnames)
  
    vals <- vector(mode="list", length=D)
    vals[[margin]] <- paste0(margin, index)
    xx <- sapply(vals, is.null)
    vals[xx] <- dimnames(X)[xx]
  
    mat <- as.matrix(expand.grid(vals))
    Y <- X[mat]
    return(array(Y,dim=p[-margin]))
  }
}

# Perform OLS on each unit of outcome separately
vox.ols.reg <- function(Y.train, X.train, find.signif=F, intercept=T,mult.adj=1) {
  N <- dim(Y.train)[1]
  p <- dim(Y.train)[-1]
  K <- ncol(X.train)
  Gamma.allk.ols <- array(NA,dim=c(K,prod(p)))
  if (find.signif) {Gamma.allk.signif <- array(NA,dim=c(K,prod(p)))}
  Y.train.vec <- array(Y.train, dim=c(N,prod(p)))
  for (vox in c(1:prod(p))) {
    y.vox <- Y.train.vec[,vox]
    if (sum(!is.na(y.vox)) > (K+1)) {
    if (intercept) {
      reg.vox <- lm(y.vox ~ X.train)
    } else {
      reg.vox <- lm(y.vox ~ 0 + X.train)
    }
    reg.sum <- summary(reg.vox)
    coefs <- reg.sum$coefficients[,1]
    signif <- reg.sum$coefficients[,4] < 0.05/mult.adj
    if (length(coefs) < K) {
      coefs <- rep(0,K)
      signif <- rep(F,K)
    }
    Gamma.allk.ols[,vox] <- coefs
    if (find.signif) {Gamma.allk.signif[,vox] <- signif}
  
    }
  }
    
  if (find.signif) {return(list(Gamma.allk.ols, Gamma.allk.signif))}
  else {return(Gamma.allk.ols)}
}

# Perform voxel-wise OLS with random intercepts (lmer)
vox.ols.reg.ri <- function(Y.train, X.train, ID, find.signif=F, mult.adj=1) {
  N <- dim(Y.train)[1]
  p <- dim(Y.train)[-1]
  K <- ncol(X.train)
  nID <- length(unique(ID))
  Gamma.allk.ols <- array(0,dim=c(K,prod(p)))
  Bi.ols <- array(0,dim=c(nID,prod(p)))
  if (find.signif) {
    Gamma.allk.signif <- array(0,dim=c(K,prod(p)))
  }
  Y.train.vec <- array(Y.train, dim=c(N,prod(p)))
  for (vox in c(1:prod(p))) {
    y.vox <- Y.train.vec[,vox]
    if (sum(!is.na(y.vox))>(K+length(unique(ID))) & sd(y.vox,na.rm=T)!=0) {
      reg.vox <- lmer(y.vox ~ X.train + (1|ID))
      if (isSingular(reg.vox)) {
        reg.vox.nori <- lm(y.vox ~ 0+X.train)
        if (length(coef(reg.vox.nori))==K) {
          Gamma.allk.ols[,vox] <- as.numeric(coef(reg.vox.nori))
        }
        if (find.signif) {
          reg.sum.nori <- summary(reg.vox.nori)
          signif <- reg.sum.nori$coefficients[,4] < 0.05 / mult.adj
          Gamma.allk.signif[,vox] <- signif
        }
      }
      if (!isSingular(reg.vox) & length(coef(reg.vox))==(K+1)) {
        coefs <- coef(reg.vox)
        Gamma.allk.ols[,vox] <- as.numeric(coefs$ID[1,2:(K+1)])
        for (i in 1:nID) {
          Bi.ols[i,vox] <- as.numeric(coefs$ID[i,1])
        }
        if (find.signif) {
          reg.sum <- summary(reg.vox)
          signif <- reg.sum$coefficients[,4] < 0.05 / mult.adj
          Gamma.allk.signif[,vox] <- signif[2:(K+1)]
        }
      }
    }
  }
  if (find.signif) {
    return(list(Gamma.allk.ols, Gamma.allk.signif, Bi.ols))
  }
  else {
    return(Gamma.allk.ols, Bi.ols)
  }
}

# Given a tensor and chosen rank, find low-rank approximation
tensor.decomp <- function(Gamma, rank) {
  require(rTensor)
  Gamma_t <- as.tensor(Gamma)
  decomp.margins <- cp(Gamma_t, num_components = rank)
  return(decomp.margins$U)
}

# given samples of tensor margins (list) returns voxel-level coefs at each mcmc iteration
getTensorCoef <- function(tensor_margins) {
  D <- length(tensor_margins)
  N_iter <- dim(tensor_margins[[1]])[1]
  p <- sapply(1:D, function(x) dim(tensor_margins[[x]])[2])
  R <- dim(tensor_margins[[1]])[3]
  K <- dim(tensor_margins[[1]])[4]
  for (d in 1:D) {
    p[d] <- dim(tensor_margins[[d]])[2]
  }
  Coefs_mcmc <- array(dim=c(N_iter,K,prod(p)))
  for (s in 1:N_iter) {
    for (k in 1:K) {
      Coefs_mcmc[s,k,] <- c(TP.rankR(lapply(tensor_margins, function(x) x[s,,,k])))
    }
  }
  return(Coefs_mcmc)
}

getR2 <- function(Y.train, X.train, ID=NA, Gamma_est, Bi_est=0) {
  N <- dim(Y.train)[1]
  nvox <- prod(dim(Y.train)[-1])
  if (all(is.na(ID))) {
    ID <- rep(1,N)
    Bi_est <- array(0,dim=c(1,nvox))
  }
  Yvec <- array(Y.train, dim=c(dim(Y.train)[1],nvox))
  
  Yhat <- t(sapply(1:N, function(x) X.train[x,] %*% Gamma_est + Bi_est[ID[x],]))
  Resid <- Yvec - Yhat
  R2_vox <- sapply(1:nvox, function(x) 1-sum(Resid[,x]^2, na.rm=T)/sum((Yvec[,x]-mean(Yvec[,x],na.rm=T))^2,na.rm=T))
  return(R2_vox)
}


