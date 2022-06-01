sample_prior <- function(prior, T_all, p_tmp, sampleDynamics = F){
  
  # linear dynamics: b, A, Q
  theta <- list()
  theta$A <- diag(p_tmp + 1)
  theta$b <- rep(0, p_tmp + 1)
  theta$Q <- diag(p_tmp + 1)*prior$sig20
  
  # mean & latent: mu, X
  muX <- matrix(Inf, nrow = p_tmp +1, ncol = T_all)
  muX[,1] <- rmvnorm(1, mean = c(prior$mu0, prior$x0),
                     sigma = diag(c(prior$Sigmu0, diag(prior$Q0))))
  
  if(sampleDynamics){
    for(k in 1:(p_tmp+1)){
      theta$Q[k,k] <- 1/rgamma(1, shape = prior$nu0/2,
                               rate = prior$nu0*prior$sig20/2)
      bAsamp <- rmvnorm(1, mean = prior$BA0,
                        sigma = theta$Q[k,k]*solve(prior$Lamb0))
      theta$b[k] <- bAsamp[1]
      theta$A[k,k] <- bAsamp[2]
    }
  }
  
  
  for(i in 1:(p_tmp+1)){
    while(sum(muX[i,] > 2) > 1){
      k <- ceiling(runif(1)*25)+10
      muX[i,] <- spline(seq(0,1, length.out = k),
                        rnorm(k)*0.3, n = T_all)$y
    }
  }
  
  
  theta$mu <- muX[1,]
  theta$X <- matrix(muX[-1,], ncol = T_all)
  
  muBar <- mean(theta$mu)
  theta$mu <- theta$mu - muBar
  
  Xbar <- apply(theta$X, 1, mean)
  theta$X <- theta$X - Xbar
  
  return(theta)
}


# matplot(t(theta$X))
