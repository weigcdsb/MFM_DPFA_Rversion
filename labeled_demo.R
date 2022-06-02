library(Rcpp)
library(rstanarm)
# library(bpr)
library(INLA)
library(pgdraw)
library(mvtnorm)
#############################################
#### load functions
# usr_dir <- 'C:\\Users\\gaw19004\\Documents\\GitHub\\MFM_DPFA_Rversion'
usr_dir <- 'D:\\github\\MFM_DPFA_Rversion'

setwd(paste0(usr_dir,'\\function'))
source('sample_prior.R')
source('update_clusParam_PG.R')

## make stand-alone Rcpp function later...
cppFunction('arma::mat FFBS(arma::mat A,
              arma::vec b,
              arma::vec m0,
              arma::mat V0,
              arma::mat Sig,
              arma::mat Yhat_tmp,
              arma::mat X_tmp,
              arma::mat w_tmp,
              int T,
              int p){
              
              arma::vec Tvec = arma::linspace(0, T-1, T);
              arma::mat m_tmp(p,T);
              arma::cube V_tmp(p,p,T);
              
              // return
              arma::mat BETA_b(p,T);
              
              
              arma::mat Kt;
              arma::mat m_tt_1;
              arma::mat V_tt_1;
              
              // (1) FF: forward-sampling
              for(int t = 0; t < T; ++t){
              
                // prediction
                if(t == 0){
                  m_tt_1 = A*m0 + b;
                  V_tt_1 = A*V0*A.t() + Sig;
                }else{
                  m_tt_1 = A*m_tmp.col(t-1) + b;
                  V_tt_1 = A*V_tmp.slice(t-1)*A.t() + Sig;
                }
              
                arma::uvec idx = arma::find_finite(Yhat_tmp.col(t));
                arma::mat X_tmp2 = X_tmp.rows(idx);
                arma::vec w_tmp2 = w_tmp.submat(idx,find(Tvec == t));
                arma::mat Yhat_tmp2 = Yhat_tmp.submat(idx,find(Tvec == t));
              
                // filtering
                Kt = V_tt_1 * X_tmp2.t()*
                  inv(X_tmp2 * V_tt_1 *X_tmp2.t() + diagmat(1/w_tmp2));
                m_tmp.col(t) = m_tt_1 + Kt*(Yhat_tmp2 - X_tmp2*m_tt_1);
                V_tmp.slice(t) = (arma::eye(p,p) - Kt*X_tmp2)*V_tt_1;
                V_tmp.slice(t) = (V_tmp.slice(t) + V_tmp.slice(t).t())/2;
              }
              
              // (2) BS: backward sampling
              BETA_b.col(T-1) = mvrnormArma(1,
              m_tmp.col(T-1),V_tmp.slice(T-1));
              
              for(int t = T-2; t >= 0; t--){
              
              arma::mat Jt = V_tmp.slice(t) *A.t()*
              inv(A*V_tmp.slice(t)*A.t() + Sig);
              arma::vec mstar_tmp = m_tmp.col(t) + 
              Jt*(BETA_b.col(t+1) - A*m_tmp.col(t) - b);
              arma::mat Vstar_tmp = (arma::eye(p,p) - Jt*A)*V_tmp.slice(t);
              Vstar_tmp = (Vstar_tmp + Vstar_tmp.t())/2;
              BETA_b.col(t) = mvrnormArma(1,mstar_tmp,Vstar_tmp);
              }
              
              
              return BETA_b;
              }',depends="RcppArmadillo",
            includes = 'arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::trans(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
}')


#############################################
# #### read data & ground truth
setwd(paste0(usr_dir,'\\data_gen'))

## observation
Y <- as.matrix(read.csv('Y.csv', header = F))
lab <- c(as.matrix(read.csv('lab.csv', header = F)))

## ground truth
delta <- c(as.matrix(read.csv('delta.csv', header = F)))
mu <- as.matrix(read.csv('mu.csv', header = F))
X <- as.matrix(read.csv('X.csv', header = F))
logLam <- as.matrix(read.csv('logLam.csv', header = F))

# for(i in 1:10){
#   Y[i,sample(1:1000, 500)] <- NA
# }

###########################################
#### pre-MCMC: priors & settings

N <- dim(Y)[1]
T_all <- dim(Y)[2]
p_max <- N - 1
nClus <- length(unique(lab))

## priors
prior <- list()
prior$mu0 <- 0
prior$Sigmu0 <- 1
prior$BA0 <- c(0,1)
prior$Lamb0 <- diag(2)
prior$nu0 <- 1
prior$sig20 <- 0.01

## pre-allocation
ng <- 1000 # total iterations

delta_fit <- matrix(NA, nrow = N, ncol = ng)
C_fit <- array(0, dim = c(N,p_max, ng))

THETA <- list()
THETA[[1]] <- list()
## initialization
delta_fit[,1] <- rnorm(N)
for(kk in 1:nClus){
  idx_tmp <- (lab == kk)
  N_tmp <- sum(idx_tmp)
  p_tmp <- 1 # 1
  
  prior$x0 <- rep(0,p_tmp)
  prior$Q0 <- diag(p_tmp)
  
  if(p_tmp > 0){
    C_fit[idx_tmp, 1:p_tmp, 1] <- matrix(rnorm(N_tmp*p_tmp), nrow = N_tmp)
  }else{
    C_fit[idx_tmp, 1, 1] <- 0
  }
  
  THETA[[1]][[kk]] <- list()
  THETA[[1]][[kk]] <- sample_prior(prior, T_all, p_tmp)
  THETA[[1]][[kk]]$p <- p_tmp
}


###########################################
#### MCMC

r_all <- 10
ACC_trace <- matrix(NA, nrow = nClus, ncol = ng)
ACC <- rep(NA, nClus)

birth_rate <- 1e-3
birth_time <- 1
alpha <- 2

for(g in 2:ng){
  
  cat(paste("iter = ", g))
  
  THETA[[g]] <- THETA[[g-1]]
  delta_fit[,g] <- delta_fit[,g-1]
  C_fit[,,g] <- C_fit[,,g-1]
  
  # (1) update cluster parameters, by (PG + FFBS + MH)
  for(j in 1:nClus){
    idx_tmp <- (lab == j)
    
    p_tmp <- THETA[[g]][[j]]$p
    prior$x0 <- rep(0,p_tmp)
    prior$Q0 <- diag(p_tmp)
    
    if(p_tmp > 0){
      C_tmp <- C_fit[idx_tmp,1:p_tmp,g-1]
    }else{
      C_tmp <- numeric(0)
    }
    
    res <- update_clusParam_PG(
      THETA[[g-1]][[j]],
      delta_fit[idx_tmp,g-1],
      C_tmp,
      Y[idx_tmp,],
      prior,
      matrix(r_all, nrow = sum(idx_tmp), ncol = T_all))
    
    THETA[[g]][[j]] <- res$theta_b
    ACC_trace[j,g] <- res$acc
    ACC[j] <- sum(ACC_trace[j,2:g])/(g-1)
  }
  
  # (2) update non-cluster parameters: intercept + loading
  for(ii in 1:N){
    lab_tmp <- lab[ii]
    Tidx_tmp <- !is.na(Y[ii,])
    
    Y_tmp <- unname(c(Y[ii,Tidx_tmp]))
    X_tmp <- t(rbind(1, THETA[[g]][[lab_tmp]]$X[,Tidx_tmp]))
    offset_tmp <- THETA[[g]][[lab_tmp]]$mu[Tidx_tmp]
    p_tmp <- THETA[[g]][[lab_tmp]]$p
    if(p_tmp > 0){
      c_tmp <- C_fit[ii,1:p_tmp,g-1]
    }else{
      c_tmp <- numeric(0)
    }
    # use stan
    deltC0 <- c(delta_fit[ii,g-1], c_tmp)
    d_tmp <- data.frame(Y = Y_tmp, X = X_tmp)
    
    fit <- stan_glm(Y ~ .-1,
                    data = d_tmp,
                    family = poisson(link="log"),
                    prior = normal(0, 1),
                    refresh = 0,
                    offset = offset_tmp,
                    chains = 1, iter = 5,
                    init = deltC0)
    deltC_samp <- c(tail(as.matrix(fit), 1))
    delta_fit[ii,g] <- deltC_samp[1]
    
    if(p_tmp > 0){
      C_fit[ii,1:p_tmp,g] <- deltC_samp[-1]
    }else{
      C_fit[ii, 1, g] <- 0
    }
  }
  
  # delta_fit[,g]
  # C_fit[,,g]
  
  
  # (3) projection
  for(j in 1:nClus){
    obsIdx <- (lab == j)
    N_tmp <- sum(obsIdx)
    p_tmp <- THETA[[g]][[j]]$p
    
    if(p_tmp > 0){
      C_tmp <- as.matrix(C_fit[obsIdx,1:p_tmp,g])
    }else{
      C_tmp <- numeric(0)
    }
    
    # for mu
    muBar <- mean(THETA[[g]][[j]]$mu)
    THETA[[g]][[j]]$mu <- THETA[[g]][[j]]$mu - muBar
    THETA[[g]][[j]]$b[1] <-
      (THETA[[g]][[j]]$A[1,1] - 1)*muBar + THETA[[g]][[j]]$b[1]
    delta_fit[obsIdx,g] <- delta_fit[obsIdx,g] - muBar
    
    # for X
    Xbar <- apply(THETA[[g]][[j]]$X,1, mean)
    THETA[[g]][[j]]$X <- THETA[[g]][[j]]$X - Xbar
    THETA[[g]][[j]]$b[-1] <- 
      (THETA[[g]][[j]]$A[-1,-1] - diag(p_tmp)) %*% Xbar +
      THETA[[g]][[j]]$b[-1]
    delta_fit[obsIdx,g] <- delta_fit[obsIdx,g] +
      c(C_tmp %*% Xbar)
  }
  
  # (4) birth-death
  for(j in 1:nClus){
    
    t_fa <- 0
    obsIdx <- (lab == j)
    N_tmp <- sum(obsIdx)
    
    while(t_fa < birth_time){
      
      p_tmp <- THETA[[g]][[j]]$p
      
      Y_bd <- c(t(Y[obsIdx,]))
      delta_bd <- delta_fit[obsIdx,g] %x% rep(1,T_all)
      mu_bd <- rep(THETA[[g]][[j]]$mu, T_all)
      
      if(p_tmp > 0){
        mll_del <- rep(0,p_tmp)
        # marginal likelihood deleting each column
        for(kk in 1:p_tmp){
          colSelect <- setdiff(1:p_tmp, kk)
          
          X_bd1 <- t(matrix(THETA[[g]][[j]]$X[colSelect,], ncol = T_all))
          X_bd <- rep(1, N_tmp) %x% X_bd1
          d_tmp <- as.data.frame(list(y = Y_bd, x = X_bd))
          
          f <- formula(paste("y ~",
                             paste(names(d_tmp)[-1],
                                   collapse='+'),
                             "-1"))
          
          fit_tmp <- inla(f,family = "poisson",
                          control.compute=list(mlik = T),
                          data = d_tmp,
                          offset = delta_bd+mu_bd,
                          control.fixed
                          =list(mean=0, prec=1))
          mll_del[kk] <- fit_tmp$mlik[1]
        }
        
        X_bd1 <- t(THETA[[g]][[j]]$X)
        X_bd <- rep(1, N_tmp) %x% X_bd1
        d_tmp <- as.data.frame(list(y = Y_bd, x = X_bd))
        f <- formula(paste("y ~",
                           paste(names(d_tmp)[-1],
                                 collapse='+'),
                           "-1"))
        
        fit_tmp <- inla(f,family = "poisson",
                        control.compute=list(mlik = T),
                        data = d_tmp,
                        offset = delta_bd+mu_bd,
                        control.fixed
                        =list(mean=0, prec=1))
        mll_all <- fit_tmp$mlik[1]
        
        delta_j <- exp(mll_del + log(birth_rate) - mll_all - log(alpha))
        delta_j[is.infinite(delta_j)] <- .Machine$double.xmax
      }else{
        delta_j <- 0
      }
      
      delta_all <- sum(delta_j)
      s <- rexp(1, rate = birth_rate + delta_all)
      t_fa <- t_fa + s
      
      if(rbinom(1,1,birth_rate/(birth_rate + delta_all)) == 1){
        if(p_tmp < p_max){
          
          # give a birth
          print('birth')
          prior$x0 <- rep(0,1)
          prior$Q0 <- diag(1)
          theta_tmp <- sample_prior(prior, T_all, 1)
          
          THETA[[g]][[j]]$b <- c(THETA[[g]][[j]]$b, theta_tmp$b[2])
          THETA[[g]][[j]]$A <- bdiag(THETA[[g]][[j]]$A, theta_tmp$A[2,2])
          THETA[[g]][[j]]$Q <- bdiag(THETA[[g]][[j]]$Q, theta_tmp$Q[2,2])
          THETA[[g]][[j]]$X <- rbind(THETA[[g]][[j]]$X, theta_tmp$X)
          
          C_fit[obsIdx, p_tmp+1, g] <- rnorm(N_tmp)
          THETA[[g]][[j]]$p <- THETA[[g]][[j]]$p + 1
        }
      }else{
        
        # give a death
        print('death')
        del_col <- c(rmultinom(1, 1, delta_j/delta_all))
        
        del_col_expand <- c(0,del_col)
        THETA[[g]][[j]]$b <- THETA[[g]][[j]]$b[!del_col_expand]
        THETA[[g]][[j]]$A <- THETA[[g]][[j]]$A[!del_col_expand,!del_col_expand]
        THETA[[g]][[j]]$Q <- THETA[[g]][[j]]$Q[!del_col_expand,!del_col_expand]
        THETA[[g]][[j]]$X <- matrix(THETA[[g]][[j]]$X[!del_col,], ncol = T_all)
        
        C_tmp <- 0*C_fit[obsIdx,,g]
        if(p_tmp > 1){
          C_tmp[,1:(p_tmp-1)] <- C_fit[obsIdx, which(!del_col), g]
        }
        
        C_fit[obsIdx,,g] <- C_tmp
        THETA[[g]][[j]]$p <- THETA[[g]][[j]]$p-1
      }
    }
    
  }
  
}





