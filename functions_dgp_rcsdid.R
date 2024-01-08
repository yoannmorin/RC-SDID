## -----------------------------------------------------------------------------
## functions_dgp_rcsdid.R
## Yoann Morin
## 
## Compute data for a DGP with repeated cross sectional data with latent factors at the group level
## -----------------------------------------------------------------------------

# Packages
library("dplyr")

#' Compute data for a DGP with repeated cross sectional data with latent factors at the group level
#' 
#' @param nco - number of control groups
#' @param ntr - number of treated groups
#' @param nT - number of periods
#' @param nT_post - number of post-treatment periods
#' @param rc_dim - list of element to simulate the repeated cross section dimension of the data. 
#' Example : rc_dim = c("scalenor", "4", "100") uses 100 baseline observations per cross sections, with a random scaling factor S 
#' drawn from a uniform law from 1 to 4, rounded to the nearest integer.
#' @param path - path to files containing FE/IFE
#' @param nbf - number of IFE. Always include group and time fixed effects.
#' @param tre - value of the treatment effect
#' @param cor_scale - correlation between individual fixed effect and the scale parameter. Either 0, 02 for 0.2, 05, 08 or 1.
#' @param w - overlap parameter
dgp_rcsdid <- function(nco, ntr, nT, nT_post, rc_dim, path, nbf, tre, cor_scale, w){
 
  # Simulate FE/IFE
    ## Import IFE/FE
  setwd(path)
  
  dfei0 = readRDS("dfei0_unif.rds")
  dfet = readRDS("dfet_norm.rds")
  dvecL0 = readRDS("dvecL0_unif.rds")
  dvecF = readRDS("dvecF_norm.rds")
  
    ## Treated group and period
  vecD = c(rep(1,ntr),rep(0,nco))
  vecT = c(rep(0,nT_post),rep(1,nT-nT_post))
  
    ## FE and IFE
  fei0 = as.matrix(c(dfei0[1:ntr,1],dfei0[(nrow(dfei0)-(nco+ntr)+ntr+1):nrow(dfei0),1]),(nco+ntr),1)
  fet = as.matrix(dfet[1:nT,],nT,1)
  
  if (nbf!=0) {
    vecL0	= as.matrix(rbind(dvecL0[1:ntr,],dvecL0[(nrow(dvecL0)-(nco+ntr)+ntr+1):nrow(dvecL0),]),(nco+ntr),ncol(dvecL0))
    vecF = as.matrix(dvecF[1:nT,],nT,nbf)
    
    fei	= as.matrix(fei0,(nco+ntr),1)
    fei[1:ntr,1] = fei[1:ntr,1] + matrix((1-w)*2*sqrt(3),ntr,1)  
    vecL = as.matrix(vecL0[,1:nbf],(nco+ntr),nbf) 
    vecL[1:ntr,] = vecL[1:ntr,] + matrix((1-w)*2*sqrt(3),ntr,nbf)
    
    ## Create dataset with IFE/FE
    dfe = data.frame(
      group = rep(1:(nco+ntr), each=length(1:nT)),
      time = rep(1:nT, length(1:(nco+ntr))),
      fei =as.matrix(fei %x% rep(1,nT),(nco+ntr)*nT,1),
      fet = as.matrix(rep(1,(nco+ntr)) %x% fet,(nco+ntr)*nT,1),
      ife = matrix(vecF[,1:nbf] %*% t(vecL),(nco+ntr)*nT,1),
      treatment = as.matrix(vecD %x% vecT))
    
    # Simulate repeated cross sectional dimension
      # Import RC random variables
    name_scale = paste("scale_param_cor",paste(cor_scale,paste(paste("w", w, sep=""),"rds", sep="."), sep="_"), sep="")
    scale_param = readRDS(name_scale)
    normevo = readRDS("normevo.rds")
    
    nscale = paste("s1", rc_dim[2], sep="_")
    nvrc = paste("n", rc_dim[3], sep="")
    
      # Create number of obs by RC
    dfe$scale = rep(scale_param[1:(nco+ntr), nscale], each=length(1:nT))
    dfe$evo = normevo[1:((nco+ntr)*(nT)), nvrc]
    
    dfe$nb0 = ifelse(dfe$time==1, round(dfe$scale)*(as.numeric(rc_dim[3])) + round(round(dfe$scale)*dfe$evo), 0)
    for (i in 1:(nco+ntr)) { 
      for (t in 2:nT) { 
        dfe$nb0[dfe$group==i & dfe$time==t] = dfe$nb0[dfe$group==i & dfe$time==(t-1)] + round(round(dfe$scale[dfe$group==i & dfe$time==(t)])*dfe$evo[dfe$group==i & dfe$time==(t)])
      }}
    
      # Transform data to RC format
    df = data.frame(
      group = rep.int(rep(1:(nco+ntr), each = length(1:nT)), dfe$nb0),
      time = rep.int(rep(1:nT, length(1:(nco+ntr))), dfe$nb0)
    )
    
    df = merge(df, dfe, by=c("group", "time"))
  }
  
  if (nbf==0) {
    fei	= as.matrix(fei0,(nco+ntr),1)
    fei[1:ntr,1] = fei[1:ntr,1] + matrix((1-w)*2*sqrt(3),ntr,1)  
    
    ## Create dataset with IFE/FE
    dfe = data.frame(
      group = rep(1:(nco+ntr), each=length(1:nT)),
      time = rep(1:nT, length(1:(nco+ntr))),
      fei =as.matrix(fei %x% rep(1,nT),(nco+ntr)*nT,1),
      fet = as.matrix(rep(1,(nco+ntr)) %x% fet,(nco+ntr)*nT,1),
      treatment = as.matrix(vecD %x% vecT))
    
    # Simulate repeated cross sectional dimension
      # Import RC random variables
    name_scale = paste("scale_param_cor",paste(cor_scale,paste(paste("w", w, sep=""),"rds", sep="."), sep="_"), sep="")
    scale_param = readRDS(name_scale)
    normevo = readRDS("normevo.rds")
    
    nscale = paste("s1", rc_dim[2], sep="_")
    nvrc = paste("n", rc_dim[3], sep="")
    
      # Create number of obs by RC
    dfe$scale = rep(scale_param[1:(nco+ntr), nscale], each=length(1:nT))
    dfe$evo = normevo[1:((nco+ntr)*(nT)), nvrc]
    
    dfe$nb0 = ifelse(dfe$time==1, round(dfe$scale)*(as.numeric(rc_dim[3])) + round(round(dfe$scale)*dfe$evo), 0)
    for (i in 1:(nco+ntr)) { 
      for (t in 2:nT) { 
        dfe$nb0[dfe$group==i & dfe$time==t] = dfe$nb0[dfe$group==i & dfe$time==(t-1)] + round(round(dfe$scale[dfe$group==i & dfe$time==(t)])*dfe$evo[dfe$group==i & dfe$time==(t)])
      }}
    
      # Transform data to RC format
    df = data.frame(
      group = rep.int(rep(1:(nco+ntr), each = length(1:nT)), dfe$nb0),
      time = rep.int(rep(1:nT, length(1:(nco+ntr))), dfe$nb0)
    )
    
    df = merge(df, dfe, by=c("group", "time"))
  }
  
  # Compute outcome and remove useless variables from the df
  if (nbf!=0) {
    df$y = as.vector(
      df$fei + df$fet + 
        tre * df$treatment + 
        df$ife)
    df$fei <- NULL
    df$fet <- NULL
    df$ife <- NULL
    df$evo <- NULL
    df$scale <- NULL
    df$nb0 <- NULL
  }
  if (nbf==0) {
    df$y = as.vector(
      df$fei + df$fet + 
        tre * df$treatment)
    df$fei <- NULL
    df$fet <- NULL
    df$evo <- NULL
    df$scale <- NULL
    df$nb0 <- NULL
  }

  # Output
  return(df)
}