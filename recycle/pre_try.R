###################### install INLA ################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)
# # 
# install.packages("INLA",
#                  repos=c(getOption("repos"),
#                          INLA="https://inla.r-inla-download.org/R/stable"),
#                  dep=TRUE)

###################### install STAN ################################
# install.packages("rstanarm")
# install.packages(c("BH", "StanHeaders", "Rcpp", "RcppEigen", "RcppParallel", "inline", "loo", "pkgbuild", "rstan"))

###################### check INLA ################################
# library(INLA)
# N <- 100
# p <- 1
# X <- matrix(rnorm(N*p), nrow = N, ncol = p)
# eta <- X*2
# y <- rpois(N, exp(eta))
# data = list(y=y, X=X)
# result <-  inla(y ~ X-1,
#                 family = "poisson",
#                 control.compute=list(mlik = T),
#                 data = data,
#                 control.fixed
#                 =list(mean=0, prec=1))
# 
# mll <- result$mlik[1]




###################### check STAN ################################
library(rstanarm)
count_data <- data.frame(
  counts = c(18,17,15,20,10,20,25,13,12),
  outcome = gl(3,1,9),
  treatment = gl(3,3)
)
initVal <- rnorm(5)

fit3 <- stan_glm(
  counts ~ outcome + treatment,
  data = count_data,
  family = poisson(link="log"),
  prior = normal(0, 2),
  refresh = 0,
  # for speed of example only
  chains = 1, iter = 2,
  init = initVal
)



trace <- as.data.frame(fit3)


for(kk in 1:10){
  fit4 <- stan_glm(
    counts ~ outcome + treatment,
    data = count_data,
    family = poisson(link="log"),
    prior = normal(0, 2),
    refresh = 0,
    # for speed of example only
    chains = 1, iter = 2,
    init = c(as.matrix(fit3))
  )
  trace <- rbind(trace, as.data.frame(fit4))
  fit3 <- fit4
}





install.packages("bpr")
library(bpr)

require(MASS) # load the data set
head(epil)

fit = sample_bpr( y ~  lbase*trt + lage + V4,
                  data = epil,
                  state = rep(0,6),
                  iter = 2)

trace2 <- fit$sim$beta[2,]


for(kk in 1:10){
  fit2 <- sample_bpr( y ~  lbase*trt + lage + V4,
                      data = epil,
                      state = fit$sim$beta[2,],
                      iter = 2, pars = list(max_dist = 1))
  
  trace2 <- rbind(trace2, fit2$sim$beta[2,])
  fit <- fit2
}


###########################################################
# install.packages('tfprobability')
# library(tfprobability)




 















