source('rcppFunctions.R')

# Functions to assign multiple values simultaneously ----------------------------

# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}



# Tensor Decomposition Functions ------------------------------------------------

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


# Assist in sampling correlation parameter alpha_dr (multi-variable case) -------------

log_g_alpha_nvars2 <- function(alpha, p, beta, w, tau, a.alpha=1, b.alpha=1,
                               D, R, K) {
  if (length(dim(beta[[1]]))<3) {
    beta <- lapply(beta, function(x) array(x,dim=c(ifelse(!is.null(dim(x)[1]),dim(x)[1],length(x)),R,K)))
  }
  if (length(dim(w))<2) {
    w <- array(w,dim=c(D,R))
  }
  log_g_alpha_dr <- array(dim=c(D,R))
  for (d in 1:D) {
    for (r in 1:R) {
      log_c_alpha = (a.alpha-1)*log(alpha[d,r]) -.5*K*(p[d]-1)*log(1-exp(-2*alpha[d,r]))
      bwb <- 0
      for (k in 1:K) {
        bwb <- bwb + (1/(tau[k]*w[d,r]*(1-exp(-2*alpha[d,r])))) * ( sum(beta[[d]][c(1,p[d]),r,k]^2)+
                                                                      (1+exp(-2*alpha[d,r]))*sum(beta[[d]][c(2:(p[d]-1)),r,k]^2)-
                                                                      2*exp(-alpha[d,r])*sum(beta[[d]][-1,r,k]*beta[[d]][-p[d],r,k]))
      }
      exp_alpha <- -.5*(bwb + 2*b.alpha*alpha[d,r])
      log_g_alpha_dr[d,r] <- log_c_alpha + exp_alpha
    }
  }
  return(log_g_alpha_dr)
}

log_g_alpha_nvars_nvisits2 <- function(alpha, p, beta, w, tau, a.alpha=1, b.alpha=1,
                                       D, R, Nt, M) {
  if (length(dim(beta[[1]]))<4) {
    beta <- lapply(beta, function(x) array(x,dim=c(ifelse(!is.null(dim(x)[1]),dim(x)[1],length(x)),R,Nt,M)))
  }
  if (length(dim(w))<2) {
    w <- array(w,dim=c(D,R))
  }
  if (length(dim(tau))<2) {
    tau <- array(tau, dim=c(Nt,M))
  }
  log_g_alpha_dr <- array(dim=c(D,R))
  for (d in 1:D) {
    for (r in 1:R) {
      log_c_alpha <- (a.alpha-1)*log(alpha[d,r]) -.5*Nt*M*(p[d]-1) * log(1-exp(-2*alpha[d,r]))
      bwb <- 0
      for (t in 1:Nt) {
        for (m in 1:M) {
          bwb <- bwb + (1/(tau[t,m]*w[d,r]*(1-exp(-2*alpha[d,r])))) * ( sum(beta[[d]][c(1,p[d]),r,t,m]^2)+
                                                                          (1+exp(-2*alpha[d,r]))*sum(beta[[d]][c(2:(p[d]-1)),r,t,m]^2)-
                                                                          2*exp(-alpha[d,r])*sum(beta[[d]][-1,r,t,m]*beta[[d]][-p[d],r,t,m]))
        }
      }
      exp_alpha <- -.5*(bwb + 2*b.alpha*alpha[d,r])
      log_g_alpha_dr[d,r] <- log_c_alpha + exp_alpha
    }
  }
  return(log_g_alpha_dr)
}



# Get coefficients/significance -------------------------------------------

# Calculates joint credible intervals using nonparametric "MDEV" method 
mdev_cred_int <- function(Gamma_k_mcmc,alpha=0.05,missing_vox=NULL) {
  if (!is.null(missing_vox)) {
    Gamma_k_mcmc[,missing_vox] <- NA
  }
  nvox <- dim(Gamma_k_mcmc)[2]
  mean_vox <- apply(Gamma_k_mcmc,2,mean)
  s_alpha_vox <- apply(Gamma_k_mcmc,2,function(x) quantile(x,c(alpha/2,1-alpha/2),na.rm = T))
  slow_alpha <- max(mean_vox-s_alpha_vox[1,],na.rm=T)
  shigh_alpha <- max(s_alpha_vox[2,]-mean_vox,na.rm=T)
  credint_vox <- rbind(mean_vox-slow_alpha,
                       mean_vox+shigh_alpha)
  return(credint_vox[1,]*credint_vox[2,]>0)
}

# Returns coefficient estimate and/or significance estimates (from credible intervals), and/or MCMC samples of coefficients
# from the MCMC samples of the tensor margins
getBTRRCoef <- function(btrr_results, term=6, burn.in=0.3, find.signif=T, signif.type="joint", alpha=0.05, median=F, output_mcmc=F, missing_vox=NULL) {
  marg <- btrr_results[[1]][[term]]
  p <- unlist(lapply(marg, function(x) dim(x)[1]))
  R <- dim(marg[[1]])[2]
  if (term != 5) {
    K <- ifelse(term%in%c(1,3),1,dim(marg[[1]])[3])
    niter <- dim(marg[[1]])[length(dim(marg[[1]]))]
    Coef_mcmc <- array(dim=c(K,length(round(burn.in*niter):niter),prod(p)))
    if (term%in%c(1,3)) {
      Coef_mcmc[1,,] <- t(sapply(round(burn.in*niter):niter, function(s) TP.rankR(lapply(marg,function(x) x[,,s]))))
    } else {
      for (k in 1:K) {
        Coef_mcmc[k,,] <- t(sapply(round(burn.in*niter):niter, function(s) TP.rankR(lapply(marg,function(x) x[,,k,s]))))
      }
    }
    if (median) {
      Coef_est <- apply(Coef_mcmc,c(1,3),median)
    } else {
      Coef_est <- apply(Coef_mcmc,c(1,3),mean)
    }
    Coef_signif <- NA*Coef_est
    if (find.signif) {
      if (signif.type=="pointwise") {
        Coef_signif <- apply(Coef_mcmc,c(1,3),function(x) quantile(x,alpha/2)*quantile(x,1-alpha/2)>0)
      } else if (signif.type=="joint") {
        for (k in 1:K) {
          Coef_signif[k,] <- mdev_cred_int(Coef_mcmc[k,,],alpha,missing_vox)
        } 
      }
    }
  } else {
    Nt <- dim(marg[[1]])[3]
    M <- dim(marg[[1]])[4]
    niter <- dim(marg[[1]])[5]
    Coef_mcmc <- array(dim=c(Nt,M,length(round(burn.in*niter):niter),prod(p)))
    for (t in 1:Nt) {
      for (m in 1:M) {
        Coef_mcmc[t,m,,] <- t(sapply(round(burn.in*niter):niter, function(s) TP.rankR(lapply(marg,function(x) x[,,t,m,s]))))
      }
    }
    if (median) {
      Coef_est <- apply(Coef_mcmc,c(1,2,4),median)
    } else {
      Coef_est <- apply(Coef_mcmc,c(1,2,4),mean)
    }
    Coef_signif <- NA*Coef_est
    if (find.signif) {
      if (signif.type=="pointwise") {
        Coef_signif <- apply(Coef_mcmc,c(1,2,4),function(x) quantile(x,alpha/2)*quantile(x,1-alpha/2)>0)
      } else if (signif.type=="joint") {
        for (t in 1:Nt) {
          for (m in 1:M) {
            Coef_signif[t,m,] <- mdev_cred_int(Coef_mcmc[t,m,,],alpha,missing_vox)
          }
        }
      }
    }
  }
  if (output_mcmc) {
    if (find.signif) {
      return(list(Coef_est, Coef_signif, Coef_mcmc))
    } else {
      return(list(Coef_est, Coef_mcmc))
    }
  } else {
    if (find.signif) {
      return(list(Coef_est, Coef_signif))
    } else {
      return(Coef_est)
    }
  }
}

# Returns the deviance of all MCMC iterations, which can be used to calculate the deviance information criteria (DIC)
# for goodness-of-fit
getDeviance_alliter <- function(btrr_results, Y.train, ID, Time, Visit, Ci, Xi, Zti, show.prog=T, prog.count=10) {
  niter <- dim(btrr_results[[5]])[2]
  terms <- !unlist(lapply(btrr_results[[1]],function(x) any(is.na(unlist(x)))))
  p <- unlist(lapply(btrr_results[[1]][[which(terms==T)[1]]],function(x) dim(x)[1]))
  Yvec <- t(apply(Y.train,1,c))
  N <- dim(Yvec)[1]
  
  if (is.null(ID)) {
    ID <- 1:N
  }
  Ni <- length(unique(ID))
  
  if (is.null(Visit)) {
    Visit <- rep(0,N)
  }
  if (is.null(Time)) {
    Time <- Visit
  }
  if (is.null(Ci)) {
    Ci <- array(0,dim=c(Ni,1))
  }
  if (is.null(Xi)) {
    Xi <- array(0,dim=c(Ni,1))
  }
  if (is.null(Zti)) {
    Zti <- array(0,dim=c(N,1))
  }
  
  Nt <- length(unique(Visit))
  M <- ncol(Ci)
  S <- ncol(Xi)
  K <- ncol(Zti)
  beta.store <- btrr_results[[1]]
  sigma2.store <- btrr_results[[5]]
  D_DIC <- rep(NA,niter)
  
  M_iter <- rep(0,prod(p))
  Bi_iter <- array(0,dim=c(Ni,prod(p)))
  Gamma_iter <- rep(0,prod(p))
  Thetai_iter <- array(0,dim=c(Ni,prod(p)))
  Btm_iter <- array(0,dim=c(Nt,M,prod(p)))
  Btm_iter_cond <- arrayC(Btm_iter,dim=c(Nt*M,prod(p)))
  Ds_iter <- array(0,dim=c(S,prod(p)))
  Gammak_iter <- array(0,dim=c(K,prod(p)))
  
  for (iter in 1:niter) {
    if (show.prog && iter %% prog.count==0) {print(iter)}
    if (terms[1]) {M_iter <- c(TP.rankR(lapply(beta.store[[1]],function(x) x[,,iter])))}
    if (terms[2]) {Bi_iter <- t(sapply(1:Ni,function(i) TP.rankR(lapply(beta.store[[2]],function(x) x[,,i,iter]))))}
    if (terms[3]) {Gamma_iter <- c(TP.rankR(lapply(beta.store[[3]],function(x) x[,,iter])))}
    if (terms[4]) {Thetai_iter <- t(sapply(1:Ni,function(i) TP.rankR(lapply(beta.store[[4]],function(x) x[,,i,iter]))))}
    if (terms[5]) {
      for (tt in 1:Nt) {
        for (mm in 1:M) {
          Btm_iter[tt,mm,] <- c(TP.rankR(lapply(beta.store[[5]],function(x) x[,,tt,mm,iter])))
        }
      }
    }
    Btm_iter_cond <- arrayC(Btm_iter,dim=c(Nt*M,prod(p)))
    if (terms[6]) {Ds_iter <- t(sapply(1:S, function(s) TP.rankR(lapply(beta.store[[6]], function(x) x[,,s,iter]))))}
    if (terms[7]) {Gammak_iter <- t(sapply(1:K,function(k) TP.rankR(lapply(beta.store[[7]],function(x) x[,,k,iter]))))}
    
    D_DIC[iter] <- getDeviance_val(Yvec,ID,Time,Visit,Ci,Xi,Zti,M_iter,Bi_iter,Gamma_iter,Thetai_iter,Btm_iter_cond,
                                   Ds_iter,Gammak_iter,sigma2.store[,iter],M)
  }
  return(D_DIC)
}


# Simulate outcome given covariates and coefficients ----------------------

# Uses full longitudinal model, but allows for lower specification (i.e. removing variables/terms)
simulate_Y <- function(p,  ID=NULL, Time=NULL, Visit=NULL, Ci=NULL, Xi=NULL, Zti=NULL, Mu=NULL, Bi=NULL, Gamma=NULL, Thetai=NULL, Btm=NULL, Ds=NULL, Gammak=NULL, noise_var=0) {
  
  # get number observations
  N <- 1
  if (!is.null(ID)) {
    N <- length(ID)
  }
  if (N==1 && !is.null(Visit)) {
    N <- length(Visit)
  }
  if (N==1 && !is.null(Time)) {
    N <- length(Time)
  }
  if (N==1 && !is.null(Zti)) {
    N <- nrow(as.matrix(Zti))
  }
  if (N==1 && !is.null(Xi)) {
    if (is.matrix(Xi)) {
      N <- nrow(Xi)
    } else if (is.vector(Xi)) {
      N <- length(Xi)
    }
  }
  
  # get number of subjects
  Ni <- N
  if (!is.null(ID)) {
    Ni <- length(unique(ID))
  }
  
  # get number of visits
  Nt <- 1
  if (!is.null(Visit)) {
    Nt <- length(unique(Visit))
  }
  
  # get covariate amounts and set defaults
  M<-S<-K<-1
  if (!is.null(Ci)) {
    M <- ncol(as.matrix(Ci))
  }
  if (!is.null(Xi)) {
    S <- ncol(as.matrix(Xi))
  }
  if (!is.null(Zti)) {
    K <- ncol(as.matrix(Zti))
  }
  
  # set default covariates and coefficients if NULL
  if (is.null(Visit)) {
    Visit <- rep(0,N)
    Time <- rep(0,N)
    Gamma<-Thetai<-Btm <- NULL
  }
  if (is.null(ID)) {
    ID <- c(1:N)
    Bi<-Thetai<-NULL
  }
  if (is.null(Ci)) {
    Ci <- array(0,dim=c(Ni,M))
    Btm <- NULL
  }
  if (is.null(Xi)) {
    Xi <- array(0,dim=c(Ni,S))
    Ds <- NULL
  }
  if (is.null(Zti)) {
    Zti <- array(0,dim=c(N,K))
    Gammak<-NULL
  }
  if (is.null(Mu)) {
    Mu <- rep(0,prod(p))
  }
  if (is.null(Bi)) {
    Bi <- array(0,dim=c(Ni,prod(p)))
  }
  if (is.null(Gamma)) {
    Gamma <- rep(0,prod(p))
  }
  if (is.null(Thetai)) {
    Thetai <- array(0,dim=c(Ni,prod(p)))
  }
  if (is.null(Btm)) {
    Btm <- array(0,dim=c(Nt,M,prod(p)))
  }
  if (is.null(Ds)) {
    Ds <- array(0,dim=c(S,prod(p)))
  }
  if (is.null(Gammak)) {
    Gammak <- array(0,dim=c(K,prod(p)))
  }
  
  # reformat covariates and coefficients
  if (min(Visit)!=0) {
    Visit <- Visit-min(Visit)
  }
  if (is.null(Time)) {
    Time <- Visit
  }
  if (Ni==1) {
    ID <- rep(1,N)
    Bi <- as.matrix(Bi)
    Thetai <- as.matrix(Thetai)
  }
  if (M==1) {
    Ci <- as.matrix(Ci)
    Btm <- array(Btm,dim=c(Nt,M,prod(p)))
  }
  if (S==1) {
    Xi <- as.matrix(Xi)
    Ds <- arrayC(Ds,dim=c(S,prod(p)))
  }
  if (K==1) {
    Zti <- as.matrix(Zti)
    Gammak <- arrayC(Gammak,dim=c(K,prod(p))) 
  }
  
  # get hat tensor
  Yhat_vec <- array(dim=c(N,prod(p)))
  for (n in 1:N) {
    Yhat_vec[n,] <- Mu + Bi[ID[n],] + Time[n]*(Gamma+Thetai[ID[n],]) +
      Ci[ID[n],]%*%Btm[Visit[n]+1,,] + Xi[ID[n],]%*%Ds + Zti[n,]%*%Gammak + rnorm(prod(p), mean=0, sd=sqrt(noise_var))
  }
  return(array(Yhat_vec,dim=c(N,p)))
}


# Single MCMC iteration for 1-variable case -------------------------------

### Performs single iteration of MCMC sampling for a single term of the longitudinal regression model
btrr.1var.1visit.1iter_rcpp <- function(Y, xc, R=2,
                                        a.sig=1, b.sig=1, a.tau=1, b.tau=1,
                                        a.lam=1, b.lam=1, a.alpha=1, b.alpha=1, sigma_log_alpha=0.01,
                                        beta.prev, w.prev, lambda.prev, tau.prev, sigma2.prev, alpha.prev,
                                        update.w=F, update.lambda=F, update.tau=T, update.sigma2=F, update.alpha=F) {
  N <- length(xc)
  p <- dim(Y)[-1]
  D <- length(p)
  Y.vec <- arrayC(Y,dim=c(N,prod(p)))
  
  beta.new <- beta.prev
  w.new <- w.prev
  lambda.new <- lambda.prev
  tau.new <- tau.prev
  sigma2.new <- sigma2.prev
  alpha.new <- alpha.prev
  
  # sample tensor margins
  if (D==3) {
    beta.new[[1]] <- sample_beta0_3D(Y.vec,xc,R,beta.prev[[1]],beta.prev[[2]],beta.prev[[3]],w.prev,tau.prev,sigma2.prev,alpha.new)
    beta.new[[2]] <- sample_beta1_3D(Y.vec,xc,R,beta.new[[1]],beta.prev[[2]],beta.prev[[3]],w.prev,tau.prev,sigma2.prev,alpha.new)
    beta.new[[3]] <- sample_beta2_3D(Y.vec,xc,R,beta.new[[1]],beta.new[[2]],beta.prev[[3]],w.prev,tau.prev,sigma2.prev,alpha.new)
  } else if (D==2) {
    beta.new[[1]] <- sample_beta0_2D(Y.vec,xc,R,beta.prev[[1]],beta.prev[[2]],w.prev,tau.prev,sigma2.prev,alpha.new)
    beta.new[[2]] <- sample_beta1_2D(Y.vec,xc,R,beta.new[[1]],beta.prev[[2]],w.prev,tau.prev,sigma2.prev,alpha.new)
  }

  # sample margin covariance (diagonal entries)
  if (update.w) {
    if (D==3) {
      c.dr_array = get_W_cdr_3D(beta.new[[1]],beta.new[[2]],beta.new[[3]],tau.prev,alpha.new)
    } else if (D==2) {
      c.dr_array = get_W_cdr_2D(beta.new[[1]],beta.new[[2]],tau.prev,alpha.new)
    }
    c.dr_array[which(c.dr_array<.Machine$double.xmin)] <- .Machine$double.xmin
    for (r in 1:R) {
      for (d in 1:D) {
        # prevent crashing due to numerical issues
        #w.new[d,r] <- max(rgig(1,1-.5*p[d],c.dr_array[d,r],lambda.prev[d,r]),.Machine$double.xmin)
        w.new[d,r] <- max(.Call("rgig",1,1-.5*p[d],c.dr_array[d,r],lambda.prev[d,r],PACKAGE = "GIGrvg"),.Machine$double.xmin)
      }
    }
  }
  
  # sample rate parameter
  if (update.lambda) {
    lambda.new <- sample_lambda(w.new,p,a.lam,b.lam)
  }
  
  # sample global variance scaling
  if (update.tau) {
    if (D==3) {
      chi_tau <- get_chiTau_3D(beta.new[[1]],beta.new[[2]],beta.new[[3]],w.new,alpha.new)
    } else if (D==2) {
      chi_tau <- get_chiTau_2D(beta.new[[1]],beta.new[[2]],w.new,alpha.new)
    }
    if (chi_tau<.Machine$double.xmin) {chi_tau<-.Machine$double.xmin}
    if (chi_tau>.Machine$double.xmax) {chi_tau<-.Machine$double.xmax}
    tau.new <- .Call("rgig",1,a.tau-.5*R*sum(p),chi_tau,2*b.tau,PACKAGE="GIGrvg")
  }
  
  # sample sigma2
  if (update.sigma2) {
    Yhat.vec <- getYhat_1var_C(xc, Gamma=TP.rankR(beta.new))
    Resid.vec <- Y.vec - Yhat.vec
    sigma2.new <- getSigma2_C(Resid.vec,a.sig,b.sig)
  }
  
  # sample alpha
  if (update.alpha) {
    alpha.cand <- arrayC(exp(rnorm(D*R, log(alpha.new), sigma_log_alpha)),dim=c(D,R))
    if (D==3) {
      acc_ratio <- exp(log_g_alpha2_C_3D(alpha.cand, p, beta.new[[1]], beta.new[[2]], beta.new[[3]], w.new, tau.new, a.alpha, b.alpha) - log_g_alpha2_C_3D(alpha.new, p, beta.new[[1]], beta.new[[2]], beta.new[[3]], w.new, tau.new, a.alpha, b.alpha))
    } else if (D==2) {
      acc_ratio <- exp(log_g_alpha2_C_2D(alpha.cand, p, beta.new[[1]], beta.new[[2]], w.new, tau.new, a.alpha, b.alpha) - log_g_alpha2_C_2D(alpha.new, p, beta.new[[1]], beta.new[[2]], w.new, tau.new, a.alpha, b.alpha))
    }
    alpha.new <- ifelse(arrayC(runif(D*R),dim=c(D,R))>acc_ratio, alpha.new, alpha.cand)
  }
  
  return(list(beta.new,w.new,lambda.new,tau.new,sigma2.new,alpha.new))
}
