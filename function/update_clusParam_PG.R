PG_FFBS <- function(BETA_a,
                    Y_tmp,X_tmp,A,b,
                    Sig,m0,
                    V0,R_tmp,
                    delta_tmp){
  
  # to debug
  # BETA_a <- rbind(theta_a$mu, theta_a$X)
  # X_tmp <- cbind(rep(1,N_tmp), C_tmp)
  # A <- theta_a$A
  # b <- theta_a$b
  # Sig <- theta_a$Q
  # m0 <- c(prior$mu0, prior$x0)
  # V0 <- diag(c(prior$Sigmu0, diag(prior$Q0)))
  
  p_tmp <- dim(BETA_a)[1]
  T_tmp <- dim(BETA_a)[2]
  N_tmp <- dim(Y_tmp)[1]
  BETA_b <- BETA_a
  
  ## (1) calculate r_{nt}, \hat{w}_{nt}, \Omega_t, \hat{y}_t
  ETA_tmp <- X_tmp %*% BETA_a + delta_tmp
  b_tmp <- R_tmp + Y_tmp
  c_tmp <- ETA_tmp - log(R_tmp)
  
  # (a) draw PG sample
  w_tmp_raw <- rep(NA, N_tmp*T_tmp)
  b_tmp_vec <- c(b_tmp)
  c_tmp_vec <- c(c_tmp)
  w_tmp_raw[!is.na(b_tmp_vec)] <-
    pgdraw(b_tmp_vec[!is.na(b_tmp_vec)],
           c_tmp_vec[!is.na(b_tmp_vec)])
  
  w_tmp <- matrix(w_tmp_raw, nrow = N_tmp)
  
  k_tmp <- (Y_tmp - R_tmp)/2 + w_tmp*
    (log(R_tmp) - delta_tmp %x% t(rep(1,T_tmp)))
  Yhat_tmp <- (1/w_tmp)*k_tmp
  
  ## (2) FF: forward filtering-- calculate m_t, V_t
  
  
  # input: 
  # (1) A, b, Sig, m0, V0
  # (2) Yhat_tmp, X_tmp, w_tmp(precision)
  # output: BETA_b
  
  BETA_b <- FFBS(A,b,m0, V0, Sig, Yhat_tmp, X_tmp,
               w_tmp, T_tmp, p_tmp)
  return(BETA_b)
}



update_clusParam_PG <- function(theta_a,
                                delta_tmp,
                                C_tmp,
                                Y_tmp,
                                prior,R_tmp){
  
  ## to debug
  # theta_a <- THETA[[g-1]][[j]]
  # delta_tmp <- delta_fit[idx_tmp,g-1]
  # C_tmp <- C_fit[idx_tmp,1:p_tmp,g-1] # C_tmp <- numeric(0)
  # Y_tmp <- Y[idx_tmp,]
  # R_tmp <- matrix(10, nrow = sum(idx_tmp), ncol = T_all)

  theta_b <- theta_a
  
  N_tmp <- dim(Y_tmp)[1]
  T_all <- dim(Y_tmp)[2]
  p_tmp <- dim(theta_a$X)[1]
  
  
  # (1) sample proposal
  muX_new <- PG_FFBS(rbind(theta_a$mu, theta_a$X),
                     Y_tmp,cbind(rep(1,N_tmp), C_tmp),
                     theta_a$A,theta_a$b,
                     theta_a$Q,c(prior$mu0, prior$x0),
                     diag(c(prior$Sigmu0, diag(prior$Q0))),
                     R_tmp,delta_tmp)
  
  # (2) MH step
  eta_ori <- cbind(rep(1, N_tmp), delta_tmp, C_tmp) %*%
    rbind(theta_a$mu, rep(1, T_all), theta_a$X)
  eta_new <- cbind(rep(1, N_tmp), delta_tmp, C_tmp) %*%
    rbind(muX_new[1,], rep(1, T_all), muX_new[-1,])
  
  lam_ori <- exp(eta_ori)
  lam_new <- exp(eta_new)
  p_ori <- 1/(1 + exp(-eta_ori + log(R_tmp)))
  p_new <- 1/(1 + exp(-eta_new + log(R_tmp)))
  
  lhr <- sum(dpois(Y_tmp, lam_new, log = T), na.rm = T) - 
    sum(dpois(Y_tmp, lam_ori, log = T), na.rm = T) +
    sum(dnbinom(Y_tmp, R_tmp, prob = p_ori, log = T), na.rm = T) - 
    sum(dnbinom(Y_tmp, R_tmp, prob = p_new, log = T), na.rm = T)
  
  if(log(runif(1)) < lhr){
    theta_b$mu <- muX_new[1,]
    theta_b$X <- matrix(muX_new[-1,], ncol = T_all)
    acc <- 1
  }else{
    theta_b$mu <- theta_a$mu
    theta_b$X <- theta_a$X
    acc <- 0
  }
  
  # (3) update linear dynamics
  muX <- rbind(theta_b$mu, theta_b$X)
  
  for(k in 1:(p_tmp+1)){
    
    tryCatch({
      # (a) update Q
      Y_BA <- muX[k,2:T_all]
      X_BA <- cbind(1, muX[k,1:(T_all-1)])
      
      Lam_n <- t(X_BA) %*% X_BA + prior$Lamb0
      BAn <- solve(Lam_n) %*% (prior$Lamb0 %*% prior$BA0 + t(X_BA) %*% Y_BA)
      
      an <- (prior$nu0 + T_all - 1)/2
      bn <- (prior$nu0*prior$sig20)/2 +
        (t(Y_BA) %*% Y_BA + 
           t(prior$BA0) %*% prior$Lamb0 %*% prior$BA0 -
           t(BAn) %*% Lam_n %*% BAn)/2
      theta_b$Q[k,k] <- 1/rgamma(1, shape = an, rate = bn)
      
      # (b) update A, b
      BAsamp <- rmvnorm(1, mean = BAn,
                        sigma = theta_b$Q[k,k]*solve(Lam_n))
      
      theta_b$b[k] <- BAsamp[1]
      theta_b$A[k,k] <- BAsamp[2]
    },
    error = function(e){})
    
  }
  return(list(theta_b = theta_b, acc = acc))
}





