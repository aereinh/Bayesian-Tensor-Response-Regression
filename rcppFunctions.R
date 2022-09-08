library(Rcpp)

# Rcpp version of R array() function
cppFunction(code='
NumericVector arrayC(NumericVector input, IntegerVector dim) {
  NumericVector output = clone(input);
  output.attr("dim") = dim;
  return output;
}
')

# Get hat tensor in 1 variable case
cppFunction(code = '
NumericMatrix getYhat_1var_C(NumericVector X, NumericVector Gamma) {
int N = X.size();
int nvox = Gamma.size();
NumericMatrix Yhat(N,nvox);
for (int n=0; n<N; ++n) {
  Yhat(n,_) = X(n)*Gamma;
}
return Yhat;
}            
')

# Sample noise variance in 1 variable case
cppFunction(code = '
NumericVector getSigma2_C(NumericMatrix Resid_vec, double a_sig, double b_sig) {
int N = Resid_vec.nrow();
NumericVector sigma2_new(N);
IntegerVector na_counts(N);
NumericVector ssr_nona(N);
NumericVector a_sig_prime(N);
NumericVector b_sig_prime(N);
for (int n=0; n<N; ++n) {
NumericVector Resid_vec_n = Resid_vec(n,_);
LogicalVector NA_indices = is_na(Resid_vec_n);
for (int i=0; i<Resid_vec_n.size();++i) {
  if (NA_indices(i)==FALSE) {
    na_counts(n)=na_counts(n)+1;
    ssr_nona(n)=ssr_nona(n)+pow(Resid_vec_n(i),2);
  }
}
a_sig_prime(n) = a_sig+.5*na_counts(n);
b_sig_prime(n) = b_sig+.5*ssr_nona(n);
sigma2_new(n) = 1/rgamma(1,a_sig_prime(n),1/b_sig_prime(n))(0);
}
return sigma2_new;
}
')

# Used to sample first dimension's tensor margin in 3D case
cppFunction(code = '
NumericMatrix sample_beta0_3D(NumericMatrix Y_vec, NumericVector xc, int R, 
                          NumericMatrix beta0_prev, NumericMatrix beta1_prev, NumericMatrix beta2_prev,
                          NumericMatrix w_prev, double tau_prev, NumericVector sigma2_prev,
                          NumericMatrix alpha_prev) {
int N = xc.size();
NumericMatrix beta0_new = clone(beta0_prev);
int p0 = beta0_prev.nrow();
int p1 = beta1_prev.nrow();
int p2 = beta2_prev.nrow();
int nvox = p0*p1*p2;
int d = 0;

// loop through r
for (int r=0; r<R; ++r) {
  // get Gamma
  NumericVector Gamma(nvox);
  for (int v=0; v<nvox; ++v) {
    int i0 = v%p0;
    int i1 = (v/p0)%p1;
    int i2 = (v/(p0*p1))%p2;
    for (int rr=0; rr<R; ++rr) {
      if (rr != r) {
        Gamma(v) = Gamma(v) + beta0_new(i0,rr)*beta1_prev(i1,rr)*beta2_prev(i2,rr);
      }
    }
  }
  
  // get Yhat_r and Yr_vec
  NumericMatrix Yhat_r(N,nvox);
  NumericMatrix Yr_vec(N,nvox);
  for (int n=0; n<N; ++n) {
    Yhat_r(n,_) = xc(n)*Gamma;
    Yr_vec(n,_) = Y_vec(n,_)-Yhat_r(n,_);
  }

  // get outer_dr (d=0)
  NumericVector outer_dr(p1*p2);
  for (int v=0; v<(p1*p2); ++v) {
    int i1 = v%p1;
    int i2 = (v/p1)%p2;
    outer_dr(v) = beta1_prev(i1,r)*beta2_prev(i2,r);
  }
  
  // loop through j
  for (int j=0; j<p0; ++j) {
  
    // get ab_drj
    double a_drj = 0;
    double b_drj = 0;
    for (int n=0; n<N; ++n) {
      NumericVector Yrvec_n = Yr_vec(n,_);
      NumericVector Y_drjn(p1*p2);
      for (int i=0; i<(p1*p2); ++i) {
        int tensor_index = j + (i%p1)*p0 + (i/p1)*p0*p1;
        Y_drjn(i) = Yrvec_n(tensor_index);
      }
      double sumSq_outer_drjn = 0;
      double sum_Y_outer_drjn = 0;
      LogicalVector NA_indices = is_na(Y_drjn);
      for (int i=0; i<Y_drjn.size(); ++i) {
        if (NA_indices(i)==FALSE) {
          sumSq_outer_drjn=sumSq_outer_drjn+pow(outer_dr(i),2);
          sum_Y_outer_drjn=sum_Y_outer_drjn+(Y_drjn(i)*outer_dr(i));
        }
      }
      a_drj = a_drj + (pow(xc(n),2)/sigma2_prev(n))*sumSq_outer_drjn;
      b_drj = b_drj + (xc(n)/sigma2_prev(n))*sum_Y_outer_drjn;
    }
    
    // get mu_drj and sigma2_drj
    double mu_drj;
    double sigma2_drj;
    if (j==0) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta0_new(j+1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else if (j==(p0-1)) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta0_new(j-1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*(beta0_new(j-1,r)+beta0_new(j+1,r)))/(a_drj*tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
    }

  //sample from normal
  beta0_new(j,r) = R::rnorm(mu_drj, sqrt(sigma2_drj));
  
  }

}

return beta0_new;

}
')

# Used to sample second dimension's tensor margin in 3D case
cppFunction(code = '
NumericMatrix sample_beta1_3D(NumericMatrix Y_vec, NumericVector xc, int R, 
                          NumericMatrix beta0_prev, NumericMatrix beta1_prev, NumericMatrix beta2_prev,
                          NumericMatrix w_prev, double tau_prev, NumericVector sigma2_prev,
                          NumericMatrix alpha_prev) {
int N = xc.size();
NumericMatrix beta1_new = clone(beta1_prev);
int p0 = beta0_prev.nrow();
int p1 = beta1_prev.nrow();
int p2 = beta2_prev.nrow();
int nvox = p0*p1*p2;
int d = 1;

// loop through r
for (int r=0; r<R; ++r) {
  // get Gamma
  NumericVector Gamma(nvox);
  for (int v=0; v<nvox; ++v) {
    int i0 = v%p0;
    int i1 = (v/p0)%p1;
    int i2 = (v/(p0*p1))%p2;
    for (int rr=0; rr<R; ++rr) {
      if (rr != r) {
        Gamma(v) = Gamma(v) + beta0_prev(i0,rr)*beta1_new(i1,rr)*beta2_prev(i2,rr);
      }
    }
  }
  
  // get Yhat_r and Yr_vec
  NumericMatrix Yhat_r(N,nvox);
  NumericMatrix Yr_vec(N,nvox);
  for (int n=0; n<N; ++n) {
    Yhat_r(n,_) = xc(n)*Gamma;
    Yr_vec(n,_) = Y_vec(n,_)-Yhat_r(n,_);
  }

  // get outer_dr (d=1)
  NumericVector outer_dr(p0*p2);
  for (int v=0; v<(p0*p2); ++v) {
    int i0 = v%p0;
    int i2 = (v/p0)%p2;
    outer_dr(v) = beta0_prev(i0,r)*beta2_prev(i2,r);
  }
  
  // loop through j
  for (int j=0; j<p1; ++j) {
  
    // get ab_drj
    double a_drj = 0;
    double b_drj = 0;
    for (int n=0; n<N; ++n) {
      NumericVector Yrvec_n = Yr_vec(n,_);
      NumericVector Y_drjn(p0*p2);
      for (int i=0; i<(p0*p2); ++i) {
        int tensor_index = j*p0 + i%p0 + (i/p0)*p0*p1;
        Y_drjn(i) = Yrvec_n(tensor_index);
      }
      double sumSq_outer_drjn = 0;
      double sum_Y_outer_drjn = 0;
      LogicalVector NA_indices = is_na(Y_drjn);
      for (int i=0; i<Y_drjn.size(); ++i) {
        if (NA_indices(i)==FALSE) {
          sumSq_outer_drjn=sumSq_outer_drjn+pow(outer_dr(i),2);
          sum_Y_outer_drjn=sum_Y_outer_drjn+(Y_drjn(i)*outer_dr(i));
        }
      }
      a_drj = a_drj + (pow(xc(n),2)/sigma2_prev(n))*sumSq_outer_drjn;
      b_drj = b_drj + (xc(n)/sigma2_prev(n))*sum_Y_outer_drjn;
    }
    
    // get mu_drj and sigma2_drj
    double mu_drj;
    double sigma2_drj;
    if (j==0) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta1_new(j+1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else if (j==(p1-1)) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta1_new(j-1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*(beta1_new(j-1,r)+beta1_new(j+1,r)))/(a_drj*tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
    }

  //sample from normal
  beta1_new(j,r) = R::rnorm(mu_drj, sqrt(sigma2_drj));
  
  }

}

return beta1_new;

}
')

# Used to sample third dimension's tensor margin in 3D case
cppFunction(code = '
NumericMatrix sample_beta2_3D(NumericMatrix Y_vec, NumericVector xc, int R, 
                          NumericMatrix beta0_prev, NumericMatrix beta1_prev, NumericMatrix beta2_prev,
                          NumericMatrix w_prev, double tau_prev, NumericVector sigma2_prev,
                          NumericMatrix alpha_prev) {
int N = xc.size();
NumericMatrix beta2_new = clone(beta2_prev);
int p0 = beta0_prev.nrow();
int p1 = beta1_prev.nrow();
int p2 = beta2_prev.nrow();
int nvox = p0*p1*p2;
int d = 2;

// loop through r
for (int r=0; r<R; ++r) {
  // get Gamma
  NumericVector Gamma(nvox);
  for (int v=0; v<nvox; ++v) {
    int i0 = v%p0;
    int i1 = (v/p0)%p1;
    int i2 = (v/(p0*p1))%p2;
    for (int rr=0; rr<R; ++rr) {
      if (rr != r) {
        Gamma(v) = Gamma(v) + beta0_prev(i0,rr)*beta1_prev(i1,rr)*beta2_new(i2,rr);
      }
    }
  }
  
  // get Yhat_r and Yr_vec
  NumericMatrix Yhat_r(N,nvox);
  NumericMatrix Yr_vec(N,nvox);
  for (int n=0; n<N; ++n) {
    Yhat_r(n,_) = xc(n)*Gamma;
    Yr_vec(n,_) = Y_vec(n,_)-Yhat_r(n,_);
  }

  // get outer_dr (d=2)
  NumericVector outer_dr(p0*p1);
  for (int v=0; v<(p0*p1); ++v) {
    int i0 = v%p0;
    int i1 = (v/p0)%p1;
    outer_dr(v) = beta0_prev(i0,r)*beta1_prev(i1,r);
  }
  
  // loop through j
  for (int j=0; j<p2; ++j) {
  
    // get ab_drj
    double a_drj = 0;
    double b_drj = 0;
    for (int n=0; n<N; ++n) {
      NumericVector Yrvec_n = Yr_vec(n,_);
      NumericVector Y_drjn(p0*p1);
      for (int i=0; i<(p0*p1); ++i) {
        int tensor_index = j*p0*p1 + (i%(p0*p1));
        Y_drjn(i) = Yrvec_n(tensor_index);
      }
      double sumSq_outer_drjn = 0;
      double sum_Y_outer_drjn = 0;
      LogicalVector NA_indices = is_na(Y_drjn);
      for (int i=0; i<Y_drjn.size(); ++i) {
        if (NA_indices(i)==FALSE) {
          sumSq_outer_drjn=sumSq_outer_drjn+pow(outer_dr(i),2);
          sum_Y_outer_drjn=sum_Y_outer_drjn+(Y_drjn(i)*outer_dr(i));
        }
      }
      a_drj = a_drj + (pow(xc(n),2)/sigma2_prev(n))*sumSq_outer_drjn;
      b_drj = b_drj + (xc(n)/sigma2_prev(n))*sum_Y_outer_drjn;
    }
    
    // get mu_drj and sigma2_drj
    double mu_drj;
    double sigma2_drj;
    if (j==0) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta2_new(j+1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else if (j==(p2-1)) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta2_new(j-1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*(beta2_new(j-1,r)+beta2_new(j+1,r)))/(a_drj*tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
    }

  //sample from normal
  beta2_new(j,r) = R::rnorm(mu_drj, sqrt(sigma2_drj));
  
  }

}

return beta2_new;
}
')

# Used to sample first dimension's tensor margin in 2D case
cppFunction(code = '
NumericMatrix sample_beta0_2D(NumericMatrix Y_vec, NumericVector xc, int R, 
                          NumericMatrix beta0_prev, NumericMatrix beta1_prev,
                          NumericMatrix w_prev, double tau_prev, NumericVector sigma2_prev,
                          NumericMatrix alpha_prev) {
int N = xc.size();
NumericMatrix beta0_new = clone(beta0_prev);
int p0 = beta0_prev.nrow();
int p1 = beta1_prev.nrow();
int nvox = p0*p1;
int d = 0;

// loop through r
for (int r=0; r<R; ++r) {
  // get Gamma
  NumericVector Gamma(nvox);
  for (int v=0; v<nvox; ++v) {
    int i0 = v%p0;
    int i1 = (v/p0)%p1;
    for (int rr=0; rr<R; ++rr) {
      if (rr != r) {
        Gamma(v) = Gamma(v) + beta0_new(i0,rr)*beta1_prev(i1,rr);
      }
    }
  }
  
  // get Yhat_r and Yr_vec
  NumericMatrix Yhat_r(N,nvox);
  NumericMatrix Yr_vec(N,nvox);
  for (int n=0; n<N; ++n) {
    Yhat_r(n,_) = xc(n)*Gamma;
    Yr_vec(n,_) = Y_vec(n,_)-Yhat_r(n,_);
  }

  // get outer_dr (d=0)
  NumericVector outer_dr = beta1_prev(_,r);
  
  // loop through j
  for (int j=0; j<p0; ++j) {
  
    // get ab_drj
    double a_drj = 0;
    double b_drj = 0;
    for (int n=0; n<N; ++n) {
      NumericVector Yrvec_n = Yr_vec(n,_);
      NumericVector Y_drjn(p1);
      for (int i=0; i<p1; ++i) {
        int tensor_index = j + i*p0;
        Y_drjn(i) = Yrvec_n(tensor_index);
      }
      double sumSq_outer_drjn = 0;
      double sum_Y_outer_drjn = 0;
      LogicalVector NA_indices = is_na(Y_drjn);
      for (int i=0; i<Y_drjn.size(); ++i) {
        if (NA_indices(i)==FALSE) {
          sumSq_outer_drjn=sumSq_outer_drjn+pow(outer_dr(i),2);
          sum_Y_outer_drjn=sum_Y_outer_drjn+(Y_drjn(i)*outer_dr(i));
        }
      }
      a_drj = a_drj + (pow(xc(n),2)/sigma2_prev(n))*sumSq_outer_drjn;
      b_drj = b_drj + (xc(n)/sigma2_prev(n))*sum_Y_outer_drjn;
    }
    
    // get mu_drj and sigma2_drj
    double mu_drj;
    double sigma2_drj;
    if (j==0) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta0_new(j+1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else if (j==(p0-1)) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta0_new(j-1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*(beta0_new(j-1,r)+beta0_new(j+1,r)))/(a_drj*tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
    }

  //sample from normal
  beta0_new(j,r) = R::rnorm(mu_drj, sqrt(sigma2_drj));
  
  }

}

return beta0_new;

}
')

# Used to sample second dimension's tensor margin in 2D case
cppFunction(code = '
NumericMatrix sample_beta1_2D(NumericMatrix Y_vec, NumericVector xc, int R, 
                          NumericMatrix beta0_prev, NumericMatrix beta1_prev,
                          NumericMatrix w_prev, double tau_prev, NumericVector sigma2_prev,
                          NumericMatrix alpha_prev) {
int N = xc.size();
NumericMatrix beta1_new = clone(beta1_prev);
int p0 = beta0_prev.nrow();
int p1 = beta1_prev.nrow();
int nvox = p0*p1;
int d = 1;

// loop through r
for (int r=0; r<R; ++r) {
  // get Gamma
  NumericVector Gamma(nvox);
  for (int v=0; v<nvox; ++v) {
    int i0 = v%p0;
    int i1 = (v/p0)%p1;
    for (int rr=0; rr<R; ++rr) {
      if (rr != r) {
        Gamma(v) = Gamma(v) + beta0_prev(i0,rr)*beta1_new(i1,rr);
      }
    }
  }
  
  // get Yhat_r and Yr_vec
  NumericMatrix Yhat_r(N,nvox);
  NumericMatrix Yr_vec(N,nvox);
  for (int n=0; n<N; ++n) {
    Yhat_r(n,_) = xc(n)*Gamma;
    Yr_vec(n,_) = Y_vec(n,_)-Yhat_r(n,_);
  }

  // get outer_dr (d=1)
  NumericVector outer_dr = beta0_prev(_,r);

  // loop through j
  for (int j=0; j<p1; ++j) {
  
    // get ab_drj
    double a_drj = 0;
    double b_drj = 0;
    for (int n=0; n<N; ++n) {
      NumericVector Yrvec_n = Yr_vec(n,_);
      NumericVector Y_drjn(p0);
      for (int i=0; i<p0; ++i) {
        int tensor_index = j*p0 + (i%p0);
        Y_drjn(i) = Yrvec_n(tensor_index);
      }
      double sumSq_outer_drjn = 0;
      double sum_Y_outer_drjn = 0;
      LogicalVector NA_indices = is_na(Y_drjn);
      for (int i=0; i<Y_drjn.size(); ++i) {
        if (NA_indices(i)==FALSE) {
          sumSq_outer_drjn=sumSq_outer_drjn+pow(outer_dr(i),2);
          sum_Y_outer_drjn=sum_Y_outer_drjn+(Y_drjn(i)*outer_dr(i));
        }
      }
      a_drj = a_drj + (pow(xc(n),2)/sigma2_prev(n))*sumSq_outer_drjn;
      b_drj = b_drj + (xc(n)/sigma2_prev(n))*sum_Y_outer_drjn;
    }
    
    // get mu_drj and sigma2_drj
    double mu_drj;
    double sigma2_drj;
    if (j==0) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta1_new(j+1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else if (j==(p1-1)) {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*beta1_new(j-1,r))/(a_drj*tau_prev*w_prev(d,r)+1);
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1);
    } else {
      mu_drj = (b_drj*tau_prev*w_prev(d,r) + exp(-alpha_prev(d,r))*(beta1_new(j-1,r)+beta1_new(j+1,r)))/(a_drj*tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
      sigma2_drj = (tau_prev*w_prev(d,r))/(a_drj+tau_prev*w_prev(d,r)+1+exp(-2*alpha_prev(d,r)));
    }

  //sample from normal
  beta1_new(j,r) = R::rnorm(mu_drj, sqrt(sigma2_drj));
  
  }

}

return beta1_new;

}
')

# Assists in sampling global variance scaling parameter (\tau) in 3D case
cppFunction(code = '
double get_chiTau_3D(NumericMatrix beta0_prev, NumericMatrix beta1_prev, NumericMatrix beta2_prev,
                  NumericMatrix w_prev,NumericMatrix alpha_prev) {
                  
int R = beta0_prev.ncol();
int D = 3;
long double chi_tau = 0;
NumericMatrix beta_d_prev;
for (int d=0; d<D; ++d) {
  if (d==0) {
    beta_d_prev = clone(beta0_prev);
  } else if (d==1) {
    beta_d_prev = clone(beta1_prev);
  } else {
    beta_d_prev = clone(beta2_prev);
  }
  for (int r=0; r<R; ++r) {
    long double sum_marg = 0;
    long double w_dr = w_prev(d,r);
      if (w_dr!=0) {
      for (int j=0; j<(beta_d_prev.nrow()); ++j) {
        if (j==0 | j==(beta_d_prev.nrow()-1)) {
          sum_marg = sum_marg + pow(beta_d_prev(j,r),2);
        }
        if (j>0 & j<(beta_d_prev.nrow()-1)) {
          sum_marg = sum_marg + (1+exp(-2*alpha_prev(d,r)))*pow(beta_d_prev(j,r),2);
        }
        if (j<(beta_d_prev.nrow()-1)) {
          sum_marg = sum_marg - 2*exp(-alpha_prev(d,r))*(beta_d_prev(j,r)*beta_d_prev(j+1,r));
        }
      }
      chi_tau = chi_tau + sum_marg/(w_dr*(1-exp(-2*alpha_prev(d,r))));
      //chi_tau = chi_tau + exp(log(sum_marg)-log_w_dr-log(1-exp(-2*alpha_prev(d,r))));
    }
  }
}
return chi_tau;
}')

# Assists in sampling global variance scaling parameter (\tau) in 2D case
cppFunction(code = '
double get_chiTau_2D(NumericMatrix beta0_prev, NumericMatrix beta1_prev,
                  NumericMatrix w_prev,NumericMatrix alpha_prev) {
                  
int R = beta0_prev.ncol();
int D = 2;
double chi_tau = 0;
NumericMatrix beta_d_prev;
for (int r=0; r<R; ++r) {
  for (int d=0; d<D; ++d) {
    if (d==0) {
      beta_d_prev = clone(beta0_prev);
    } else {
      beta_d_prev = clone(beta1_prev);
    }
    
    double sum_marg = 0;
    for (int j=0; j<(beta_d_prev.nrow()); ++j) {
      if (j==0 | j==(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + pow(beta_d_prev(j,r),2);
      }
      if (j>0 & j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + (1+exp(-2*alpha_prev(d,r)))*pow(beta_d_prev(j,r),2);
      }
      if (j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg - 2*exp(-alpha_prev(d,r))*(beta_d_prev(j,r)*beta_d_prev(j+1,r));
      }
    }
    
    chi_tau = chi_tau + sum_marg/(w_prev(d,r)*(1-exp(-2*alpha_prev(d,r))));
  }
}
return chi_tau;
}')

# Assists in sampling dimension and rank-specific diagonal (co)variance terms (w_dr) in 3D case
cppFunction(code='
NumericMatrix get_W_cdr_3D(NumericMatrix beta0_prev, NumericMatrix beta1_prev, NumericMatrix beta2_prev,
                  long double tau_prev,NumericMatrix alpha_prev) {
int R = beta0_prev.ncol();
int D = 3;
NumericMatrix W_cdr(D,R);
NumericMatrix beta_d_prev;
for (int r=0; r<R; ++r) {
  for (int d=0; d<D; ++d) {
    if (d==0) {
      beta_d_prev = clone(beta0_prev);
    } else if (d==1) {
      beta_d_prev = clone(beta1_prev);
    } else {
      beta_d_prev = clone(beta2_prev);
    }
    
    long double sum_marg = 0;
    for (int j=0; j<(beta_d_prev.nrow()); ++j) {
      if (j==0 | j==(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + pow(beta_d_prev(j,r),2);
      }
      if (j>0 & j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + (1+exp(-2*alpha_prev(d,r)))*pow(beta_d_prev(j,r),2);
      }
      if (j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg - 2*exp(-alpha_prev(d,r))*(beta_d_prev(j,r)*beta_d_prev(j+1,r));
      }
    }
    
    W_cdr(d,r) = sum_marg/tau_prev;
  }
}
return W_cdr;
}            
')

# Assists in sampling dimension and rank-specific diagonal (co)variance terms (w_dr) in 2D case
cppFunction(code='
NumericMatrix get_W_cdr_2D(NumericMatrix beta0_prev, NumericMatrix beta1_prev,
                  double tau_prev,NumericMatrix alpha_prev) {
int R = beta0_prev.ncol();
int D = 2;
NumericMatrix W_cdr(D,R);
NumericMatrix beta_d_prev;
for (int r=0; r<R; ++r) {
  for (int d=0; d<D; ++d) {
    if (d==0) {
      beta_d_prev = clone(beta0_prev);
    } else {
      beta_d_prev = clone(beta1_prev);
    }
    
    double sum_marg = 0;
    for (int j=0; j<(beta_d_prev.nrow()); ++j) {
      if (j==0 | j==(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + pow(beta_d_prev(j,r),2);
      }
      if (j>0 & j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + (1+exp(-2*alpha_prev(d,r)))*pow(beta_d_prev(j,r),2);
      }
      if (j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg - 2*exp(-alpha_prev(d,r))*(beta_d_prev(j,r)*beta_d_prev(j+1,r));
      }
    }
    
    W_cdr(d,r) = sum_marg/tau_prev;
  }
}
return W_cdr;
}            
')

# Samples rate parameter corresponding to diagonal (co)variance terms (\lambda_dr)
cppFunction(code = '
NumericMatrix sample_lambda(NumericMatrix w_prev, IntegerVector p, double a_lam, double b_lam) {
int R = w_prev.ncol();
int D = w_prev.nrow();
NumericMatrix lambda_new(D,R);
for (int r=0; r<R; ++r) {
  for (int d=0; d<D; ++d) {
    lambda_new(d,r) = R::rgamma(a_lam+p(d), 1/(b_lam+.5*p(d)*w_prev(d,r)));
  }
}
return lambda_new;
}            
')

# Assists in sampling correlation term (\alpha_dr) in 3D case
cppFunction(code = '
NumericMatrix log_g_alpha2_C_3D(NumericMatrix alpha_prev, IntegerVector p, 
                             NumericMatrix beta0_prev, NumericMatrix beta1_prev, NumericMatrix beta2_prev,
                             NumericMatrix w_prev, double tau_prev, double a_alpha, double b_alpha) {
int D = alpha_prev.nrow();
int R = alpha_prev.ncol();
NumericMatrix log_g_alpha(D,R);
NumericMatrix beta_d_prev;
for (int r=0; r<R; ++r) {
  for (int d=0; d<D; ++d) {
    if (d==0) {
      beta_d_prev = clone(beta0_prev);
    } else if (d==1) {
      beta_d_prev = clone(beta1_prev);
    } else {
      beta_d_prev = clone(beta2_prev);
    }
    
    double sum_marg = 0;
    for (int j=0; j<(beta_d_prev.nrow()); ++j) {
      if (j==0 | j==(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + pow(beta_d_prev(j,r),2);
      }
      if (j>0 & j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + (1+exp(-2*alpha_prev(d,r)))*pow(beta_d_prev(j,r),2);
      }
      if (j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg - 2*exp(-alpha_prev(d,r))*(beta_d_prev(j,r)*beta_d_prev(j+1,r));
      }
    }
    double bwb = sum_marg/(tau_prev*w_prev(d,r)*(1-exp(-2*alpha_prev(d,r))));
    double exp_alpha = -.5*(bwb+2*b_alpha*alpha_prev(d,r));
    double log_c_alpha = (a_alpha-1)*log(alpha_prev(d,r))-.5*(p(d)-1)*log(1-exp(-2*alpha_prev(d,r)));
    log_g_alpha(d,r) = log_c_alpha+exp_alpha;
  }
}
return log_g_alpha;
}            
')

# Assists in sampling correlation term (\alpha_dr) in 2D case
cppFunction(code = '
NumericMatrix log_g_alpha2_C_2D(NumericMatrix alpha_prev, IntegerVector p, 
                             NumericMatrix beta0_prev, NumericMatrix beta1_prev,
                             NumericMatrix w_prev, double tau_prev, double a_alpha, double b_alpha) {
int D = alpha_prev.nrow();
int R = alpha_prev.ncol();
NumericMatrix log_g_alpha(D,R);
NumericMatrix beta_d_prev;
for (int r=0; r<R; ++r) {
  for (int d=0; d<D; ++d) {
    if (d==0) {
      beta_d_prev = clone(beta0_prev);
    } else {
      beta_d_prev = clone(beta1_prev);
    }
    
    double sum_marg = 0;
    for (int j=0; j<(beta_d_prev.nrow()); ++j) {
      if (j==0 | j==(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + pow(beta_d_prev(j,r),2);
      }
      if (j>0 & j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg + (1+exp(-2*alpha_prev(d,r)))*pow(beta_d_prev(j,r),2);
      }
      if (j<(beta_d_prev.nrow()-1)) {
        sum_marg = sum_marg - 2*exp(-alpha_prev(d,r))*(beta_d_prev(j,r)*beta_d_prev(j+1,r));
      }
    }
    double bwb = sum_marg/(tau_prev*w_prev(d,r)*(1-exp(-2*alpha_prev(d,r))));
    double exp_alpha = -.5*(bwb+2*b_alpha*alpha_prev(d,r));
    double log_c_alpha = (a_alpha-1)*log(alpha_prev(d,r))-.5*(p(d)-1)*log(1-exp(-2*alpha_prev(d,r)));
    log_g_alpha(d,r) = log_c_alpha+exp_alpha;
  }
}
return log_g_alpha;
}            
')

# Gets vectorized residual of the full longitudinal model
cppFunction(code = '
NumericMatrix getResid_fullModel(NumericMatrix Yvec, IntegerVector ID, NumericVector Time, IntegerVector Visit,
                                 NumericVector Ci, NumericMatrix Xi, NumericMatrix Zti,
                                 NumericVector M_iter, NumericMatrix Bi_iter, NumericVector Gamma_iter,
                                 NumericMatrix Thetai_iter, NumericMatrix Btm_iter, NumericMatrix Ds_iter, NumericMatrix Gammak_iter,
                                 int M) {
NumericMatrix Resid_vec = clone(Yvec);
int N = Yvec.nrow();
int S = Xi.ncol();
int K = Zti.ncol();
int Nt = Btm_iter.nrow()/M;
for (int n=0; n<N; ++n) {
  Resid_vec(n,_) = Resid_vec(n,_)-M_iter-Bi_iter((ID(n)-1),_)-(Gamma_iter+Thetai_iter((ID(n)-1),_))*Time(n);
  for (int m=0; m<M; ++m) {
    Resid_vec(n,_) = Resid_vec(n,_)-Btm_iter(Visit(n)+m*Nt,_)*Ci(ID(n)-1,m);
  }
  for (int s=0; s<S; ++s) {
    Resid_vec(n,_) = Resid_vec(n,_)-Ds_iter(s,_)*Xi((ID(n)-1),s);
  }
  for (int k=0; k<K; ++k) {
    Resid_vec(n,_) = Resid_vec(n,_)-Gammak_iter(k,_)*Zti(n,k);
  }
}
return Resid_vec;
}            
')

# Returns deviance for a given set of outcomes, covariates, coefficients, and noise variances
cppFunction(code = '
double getDeviance_val(NumericMatrix Yvec, IntegerVector ID, NumericVector Time, IntegerVector Visit, 
                  NumericVector Ci, NumericMatrix Xi, NumericMatrix Zti, 
                  NumericVector M_iter, NumericMatrix Bi_iter, NumericVector Gamma_iter, 
                  NumericMatrix Thetai_iter, NumericMatrix Btm_iter,NumericMatrix Ds_iter, NumericMatrix Gammak_iter, 
                  NumericVector sigma2_iter, int M) {
                  
double dic_tot = 0;
NumericMatrix Resid_vec = clone(Yvec);
int N = Yvec.nrow();
NumericVector ssr_nona(N);
int S = Xi.ncol();
int K = Zti.ncol();
int Nt = Btm_iter.nrow()/M;
int V = Yvec.ncol();
LogicalVector NA_indices(V);
for (int n=0; n<N; ++n) {
  NA_indices = is_na(Yvec(n,_));
  Resid_vec(n,_) = Resid_vec(n,_)-M_iter-Bi_iter((ID(n)-1),_)-(Gamma_iter+Thetai_iter((ID(n)-1),_))*Time(n);
  for (int m=0; m<M; ++m) {
    Resid_vec(n,_) = Resid_vec(n,_)-Btm_iter(Visit(n)+m*Nt,_)*Ci(ID(n)-1,m);
  }
  for (int s=0; s<S; ++s) {
    Resid_vec(n,_) = Resid_vec(n,_)-Ds_iter(s,_)*Xi((ID(n)-1),s);
  }
  for (int k=0; k<K; ++k) {
    Resid_vec(n,_) = Resid_vec(n,_)-Gammak_iter(k,_)*Zti(n,k);
  }
  for (int v=0; v<V; ++v) {
    if (NA_indices(v)==FALSE) {
      ssr_nona(n) = ssr_nona(n) + pow(Resid_vec(n,v),2);
    }
  }
  dic_tot = dic_tot + (log(sigma2_iter(n))/sigma2_iter(n))*ssr_nona(n);
}
return dic_tot;
}
')

# Samples noise variance under the full longitudinal model
cppFunction(code = '
NumericVector sampleSigma_fullModel(NumericMatrix Yvec, IntegerVector ID, NumericVector Time, IntegerVector Visit,
                                 NumericVector Ci, NumericMatrix Xi, NumericMatrix Zti,
                                 NumericVector M_iter, NumericMatrix Bi_iter, NumericVector Gamma_iter,
                                 NumericMatrix Thetai_iter, NumericMatrix Btm_iter, NumericMatrix Ds_iter, NumericMatrix Gammak_iter,
                                 int M, double a_sig, double b_sig) {
NumericMatrix Resid_vec = clone(Yvec);
int N = Yvec.nrow();
int V = Yvec.ncol();
int S = Xi.ncol();
int K = Zti.ncol();
int Nt = Btm_iter.nrow()/M;
NumericVector sigma2_new(N);
for (int n=0; n<N; ++n) {
  Resid_vec(n,_) = Resid_vec(n,_)-M_iter-Bi_iter((ID(n)-1),_)-(Gamma_iter+Thetai_iter((ID(n)-1),_))*Time(n);
  for (int m=0; m<M; ++m) {
    Resid_vec(n,_) = Resid_vec(n,_)-Btm_iter(Visit(n)+m*Nt,_)*Ci(ID(n)-1,m);
  }
  for (int s=0; s<S; ++s) {
    Resid_vec(n,_) = Resid_vec(n,_)-Ds_iter(s,_)*Xi((ID(n)-1),s);
  }
  for (int k=0; k<K; ++k) {
    Resid_vec(n,_) = Resid_vec(n,_)-Gammak_iter(k,_)*Zti(n,k);
  }
  
  int notNA_n = 0;
  double SSR_n = 0;
  LogicalVector NA_indices = is_na(Resid_vec(n,_));
  for (int v=0; v<V; ++v) {
    if (NA_indices(v)==FALSE) {
      notNA_n = notNA_n+1;
      SSR_n = SSR_n + pow(Resid_vec(n,v),2);
    }
  }
  sigma2_new(n) = 1/R::rgamma(a_sig+.5*notNA_n, 1/(b_sig+.5*SSR_n));
}
return sigma2_new;
}            
')