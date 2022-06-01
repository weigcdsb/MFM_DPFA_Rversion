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
  
  Omega_tmp <- array(0,dim = c(N_tmp, N_tmp, T_tmp))
  Yhat_tmp <- matrix(0, nrow = N_tmp, ncol = T_tmp)
  
  for(t in 1:T_tmp){
    # (b) calculate \Omega_t
    Omega_tmp[,,t] <- diag(w_tmp[,t])
    
    # (c) calculate yhat_t
    kt_tmp <- (Y_tmp[,t] - R_tmp[,t])/2 +
      w_tmp[,t]*(log(R_tmp[,t]) - delta_tmp)
    Yhat_tmp[,t] <- (1/w_tmp[,t])*kt_tmp
  }
  
  ## (2) FF: forward filtering-- calculate m_t, V_t
  m_tmp <- matrix(0, nrow = p_tmp, ncol = T_tmp)
  V_tmp <- array(0, dim = c(p_tmp, p_tmp, T_tmp))
  
  for(t in 1:T_tmp){
    
    if(t == 1){
      m_tt_1 <- A %*% m0 + b
      V_tt_1 <- A %*% V0 %*% t(A) + Sig
    }else{
      m_tt_1 <-  A %*% m_tmp[,t-1] + b
      V_tt_1 <- A %*% V_tmp[,,t-1]%*% t(A) + Sig
    }
    
    obsIdx <- !is.na(Yhat_tmp[,t])
    X_tmp2 <- X_tmp[obsIdx,]
    Omega_tmp2 <- Omega_tmp[obsIdx,obsIdx,t]
    Yhat_tmp2 <- Yhat_tmp[obsIdx,t]
    
    
    if(sum(obsIdx) > 0){
      Kt <- V_tt_1 %*% t(X_tmp2) %*%
        solve(X_tmp2 %*% V_tt_1 %*% t(X_tmp2) + solve(Omega_tmp2))
      m_tmp[,t] <- as.vector(m_tt_1 + Kt %*% (Yhat_tmp2 - X_tmp2 %*% m_tt_1))
      V_tmp[,,t] <- as.matrix((diag(p_tmp) - Kt %*% X_tmp2) %*% V_tt_1)
      V_tmp[,,t] <- (V_tmp[,,t] + t(V_tmp[,,t]))/2
    }else{
      m_tmp[,t] <- m_tt_1
      V_tmp[,,t] <- V_tt_1
    }
    
  }
  
  ## (3) BS: backward sampling
  BETA_b[,T_tmp] <- rmvnorm(1, mean = m_tmp[,T_tmp],
                            sigma = as.matrix(V_tmp[,,T_tmp]))
  
  for(t in (T_tmp-1):1){
    # print(t)
    Jt <- V_tmp[,,t] %*% t(A) %*%
      solve(t(A) %*% V_tmp[,,t] %*% t(A) + Sig)
    mstar_tmp <- m_tmp[,t] + Jt %*% 
      (BETA_b[,t+1] - A %*% m_tmp[,t] - b)
    Vstar_tmp <- (diag(p_tmp) - Jt %*% A) %*% V_tmp[,,t]
    Vstar_tmp <- as.matrix((Vstar_tmp + t(Vstar_tmp))/2)
    BETA_b[,t] <- rmvnorm(1, mean = mstar_tmp,
                          sigma = Vstar_tmp)
  }
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





