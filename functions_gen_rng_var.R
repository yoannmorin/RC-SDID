## -----------------------------------------------------------------------------
## functions_gen_rng_var.R
## Yoann Morin
## 
## This creates R functions to generate random variables used in simulations
## -----------------------------------------------------------------------------

# Packages
library("rngtools")

#' Compute factors, factor loadings, individual and time fixed effects, and covariates.
#' Also generates the scaling parameter and evolution parameter to compute the number of observations in each cross section.
#' 
#' @param path - path to the location of used to store computed random variables.
functions_gen_rng_var <- function(path){
  
  setwd(path)
  
  ##################
  ##### FE/IFE #####
  ##################
  ##### Uniform
  # Create
  size = 1000
  dfei0_unif_temp = runif(n = size, min = 0, max = 1 )
  dfei0_unif = data.frame(sqrt(3) + dfei0_unif_temp*(sqrt(3) + sqrt(3)))
  dvecL0_unif =  data.frame(matrix(runif(n = size*8, min = -sqrt(3), max = sqrt(3)), size, 8))
  
  # Export
  saveRDS(dfei0_unif, file="dfei0_unif.rds")
  saveRDS(dvecL0_unif, file="dvecL0_unif.rds")
  
  
  ##### Normal
  # Create
  size = 1000
  dfet_norm = data.frame(rnorm(size,1))
  dvecF_norm = data.frame(matrix(rnorm(size*8,1), size, 8))
  
  # Export
  saveRDS(dfet_norm, file="dfet_norm.rds")
  saveRDS(dvecF_norm, file="dvecF_norm.rds")
  
  
  
  
  ##############
  ##### RC #####
  ##############
  ##### Function to correlate FE with scaling parameter
  rbvunif <- function(unif_e,rho,size) {
    old <- RNGseed() # backup old seed
    set.seed(987) # set new seed
    
    x <- unif_e
    if ((rho > 1.0) || (rho < -1.0)) {
      stop('rbvunif::rho not in [-1,+1]')
    }
    else if (rho==1.0) y <- x
    else if (rho==-1.0) y <- (1-x)
    else if (rho==0.0) y <- runif(size)
    else {
      a <- (sqrt((49+rho)/(1+rho))-5)/2
      u <- rbeta(size,a,1.0)
      y <- runif(size)
      y <- ifelse(y<0.5 ,abs(u-x), 1-abs(1-u-x) )
    }
    RNGseed(old) # restore old seed
    return(y)
  }
  
  
  ##### Scaling parameter
  function_scale_param <- function(rho, w, ntr) {
    scale_param_cor = data.frame(
      s1_1 = 1,
      s1_2 = ceiling(0 + rbvunif(dfei0_unif_temp, rho, size)*(2-0)),
      s1_4 = ceiling(0 + rbvunif(dfei0_unif_temp, rho, size)*(4-0)),
      s1_6 = ceiling(0 + rbvunif(dfei0_unif_temp, rho, size)*(6-0)),
      s1_8 = ceiling(0 + rbvunif(dfei0_unif_temp, rho, size)*(8-0)),
      s1_10 = ceiling(0 + rbvunif(dfei0_unif_temp, rho, size)*(10-0)),
      s1_15 = ceiling(0 + rbvunif(dfei0_unif_temp, rho, size)*(15-0)),
      s1_20 = ceiling(0 + rbvunif(dfei0_unif_temp, rho, size)*(20-0))
    )
    for (i in 2:(ncol(scale_param_cor))) {
      scale_param_cor[1:ntr,i] = scale_param_cor[1:ntr,i] + matrix((1-w)*2*sqrt(3),ntr,1)
      max_e = as.numeric(substr(names(scale_param_cor)[i], 4, 5))
      scale_param_cor[1:ntr,i] = ifelse(scale_param_cor[1:ntr,i]>max_e, max_e, ceiling(scale_param_cor[1:ntr,i]))
    }
    return(scale_param_cor)
  }
  
  scale_param_cor0_w0 = function_scale_param(0, 0, 1)
  scale_param_cor02_w0 = function_scale_param(0.2, 0, 1)
  scale_param_cor05_w0 = function_scale_param(0.5, 0, 1)
  scale_param_cor08_w0 = function_scale_param(0.8, 0, 1)
  scale_param_cor1_w0 = function_scale_param(1, 0, 1)
  
  scale_param_cor0_w0.2 = function_scale_param(0, 0.2, 1)
  scale_param_cor02_w0.2 = function_scale_param(0.2, 0.2, 1)
  scale_param_cor05_w0.2 = function_scale_param(0.5, 0.2, 1)
  scale_param_cor08_w0.2 = function_scale_param(0.8, 0.2, 1)
  scale_param_cor1_w0.2 = function_scale_param(1, 0.2, 1)
  
  scale_param_cor0_w0.4 = function_scale_param(0, 0.4, 1)
  scale_param_cor02_w0.4 = function_scale_param(0.2, 0.4, 1)
  scale_param_cor05_w0.4 = function_scale_param(0.5, 0.4, 1)
  scale_param_cor08_w0.4 = function_scale_param(0.8, 0.4, 1)
  scale_param_cor1_w0.4 = function_scale_param(1, 0.4, 1)
  
  scale_param_cor0_w0.6 = function_scale_param(0, 0.6, 1)
  scale_param_cor02_w0.6 = function_scale_param(0.2, 0.6, 1)
  scale_param_cor05_w0.6 = function_scale_param(0.5, 0.6, 1)
  scale_param_cor08_w0.6 = function_scale_param(0.8, 0.6, 1)
  scale_param_cor1_w0.6 = function_scale_param(1, 0.6, 1)
  
  scale_param_cor0_w0.8 = function_scale_param(0, 0.8, 1)
  scale_param_cor02_w0.8 = function_scale_param(0.2, 0.8, 1)
  scale_param_cor05_w0.8 = function_scale_param(0.5, 0.8, 1)
  scale_param_cor08_w0.8 = function_scale_param(0.8, 0.8, 1)
  scale_param_cor1_w0.8 = function_scale_param(1, 0.8, 1)
  
  scale_param_cor0_w1 = function_scale_param(0, 1, 1)
  scale_param_cor02_w1 = function_scale_param(0.2, 1, 1)
  scale_param_cor05_w1 = function_scale_param(0.5, 1, 1)
  scale_param_cor08_w1 = function_scale_param(0.8, 1, 1)
  scale_param_cor1_w1 = function_scale_param(1, 1, 1)
  
  # Export 
  saveRDS(scale_param_cor0_w0, file="scale_param_cor0_w0.rds")
  saveRDS(scale_param_cor02_w0, file="scale_param_cor02_w0.rds")
  saveRDS(scale_param_cor05_w0, file="scale_param_cor05_w0.rds")
  saveRDS(scale_param_cor08_w0, file="scale_param_cor08_w0.rds")
  saveRDS(scale_param_cor1_w0, file="scale_param_cor1_w0.rds")
  
  saveRDS(scale_param_cor0_w0.2, file="scale_param_cor0_w0.2.rds")
  saveRDS(scale_param_cor02_w0.2, file="scale_param_cor02_w0.2.rds")
  saveRDS(scale_param_cor05_w0.2, file="scale_param_cor05_w0.2.rds")
  saveRDS(scale_param_cor08_w0.2, file="scale_param_cor08_w0.2.rds")
  saveRDS(scale_param_cor1_w0.2, file="scale_param_cor1_w0.2.rds")
  
  saveRDS(scale_param_cor0_w0.4, file="scale_param_cor0_w0.4.rds")
  saveRDS(scale_param_cor02_w0.4, file="scale_param_cor02_w0.4.rds")
  saveRDS(scale_param_cor05_w0.4, file="scale_param_cor05_w0.4.rds")
  saveRDS(scale_param_cor08_w0.4, file="scale_param_cor08_w0.4.rds")
  saveRDS(scale_param_cor1_w0.4, file="scale_param_cor1_w0.4.rds")
  
  saveRDS(scale_param_cor0_w0.6, file="scale_param_cor0_w0.6.rds")
  saveRDS(scale_param_cor02_w0.6, file="scale_param_cor02_w0.6.rds")
  saveRDS(scale_param_cor05_w0.6, file="scale_param_cor05_w0.6.rds")
  saveRDS(scale_param_cor08_w0.6, file="scale_param_cor08_w0.6.rds")
  saveRDS(scale_param_cor1_w0.6, file="scale_param_cor1_w0.6.rds")
  
  saveRDS(scale_param_cor0_w0.8, file="scale_param_cor0_w0.8.rds")
  saveRDS(scale_param_cor02_w0.8, file="scale_param_cor02_w0.8.rds")
  saveRDS(scale_param_cor05_w0.8, file="scale_param_cor05_w0.8.rds")
  saveRDS(scale_param_cor08_w0.8, file="scale_param_cor08_w0.8.rds")
  saveRDS(scale_param_cor1_w0.8, file="scale_param_cor1_w0.8.rds")
  
  saveRDS(scale_param_cor0_w1, file="scale_param_cor0_w1.rds")
  saveRDS(scale_param_cor02_w1, file="scale_param_cor02_w1.rds")
  saveRDS(scale_param_cor05_w1, file="scale_param_cor05_w1.rds")
  saveRDS(scale_param_cor08_w1, file="scale_param_cor08_w1.rds")
  saveRDS(scale_param_cor1_w1, file="scale_param_cor1_w1.rds")
  
  
  ##### Evolution of obs
  # Create
  sev=0.02
  old <- RNGseed()
  set.seed(123*987)
  normevo = data.frame(n100 = rnorm(2500, mean=100*sev, sd=sqrt(100)/2))
  set.seed(123*987)
  normevo$n75 = rnorm(2500, mean=75*sev, sd=sqrt(75)/2)
  set.seed(123*987)
  normevo$n50 = rnorm(2500, mean=50*sev, sd=sqrt(50)/2)
  set.seed(123*987)
  normevo$n25 = rnorm(2500, mean=25*sev, sd=sqrt(25)/2)
  RNGseed(old)
  
  # Export
  saveRDS(normevo, file="normevo.rds")
}