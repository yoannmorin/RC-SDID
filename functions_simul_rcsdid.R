## -----------------------------------------------------------------------------
## functions_simul_rcsdid.R
## Yoann Morin
## 
## Compute simulations
## -----------------------------------------------------------------------------

# Packages
source("functions_estimate_rcsdid.R")
source("functions_dgp_rcsdid.R")
library("foreach")
library("doParallel")
library("doRNG")

#' Compute simulations for the DGP and estimate the model using different estimates.
#' @param nrep - number of repetitions for each case
#' @param vnbf - vector with each case of number of IFE to test
#' @param samp_size - a list of lists. the number of control groups (always 1 treated) in 1st position, period in 2nd position (post-treatment = period/2 so use a multiple of 2). and then list of element to simulate 
#' the repeated cross section dimension of the data. The 6th parameter corresponds to the cor_scale parameter from the dgp_rcsdid function.
#' The 7th parameter is the 
#' #' Example : samp_size = list(list("20", "20", "scalenor", "4", "100", "08", "0.6")) uses 100 baseline observations per cross sections, with a random scaling factor S 
#' drawn from a uniform law from 1 to 4, rounded to the nearest integer. The simulated data has 20 control groups and 20 periods. 
#' The correlation between the individual fixed effect and the scale parameter is 0.8 and w = 0.6 is the overlap parameter.
#' @param path - path to files containing FE/IFE
#' @param tre - value of the treatment effect
#' @param estim - "rcsdid" to compute rcsdid estimate (the default), "all" to also compute DID, sdid on aggregated data and sdid without RC adjustment.
#' @param outfile - name of the file to output simulation results, must end with ".rds"
functions_simul_rcsdid <- function(nrep, vnbf, samp_size, path, tre, estim, outfile){
  # Output results
  outres=list()
  
  for (i.ncase in 1:length(samp_size)) {
    # Creation rc_dim
    rc_dim = samp_size[[i.ncase]][ c(3, 4, 5)]
    cor_scale = samp_size[[i.ncase]][c(6)]
    w=as.numeric(samp_size[[i.ncase]][c(7)])
    
    # creation period and group size
    nco = as.numeric(samp_size[[i.ncase]][1])
    ntr = 1
    nT = as.numeric(samp_size[[i.ncase]][2])
    nT_post = round(as.numeric(samp_size[[i.ncase]][2])/2)
    
    # other GDP setup
    path = path
    tre = tre
    
    # Loop over number of IFE
    for (i.nnbf in 1:length(vnbf)) {
      nbf = vnbf[i.nnbf]
      
      # Create simulated dataset
      df = dgp_rcsdid(nco=nco, ntr=ntr, nT=nT, nT_post=nT_post, rc_dim=rc_dim, path=path, nbf=nbf, tre=tre, cor_scale=cor_scale, w=w)
      
      # Print case running
      print(paste(nco, nT, paste(rc_dim ,collapse="_"), "IFE", nbf, "cor_scale", cor_scale, "w", w, sep="_"))
      
      # Loop for simulations 
        # Parallel clusters
      cores=detectCores()
      cl <- makeCluster(16) 
      registerDoParallel(cl)
      
      res_sim <- foreach(i=1:nrep, .packages=c("fixest", "dplyr", "synthdid"), .combine=rbind) %dorng% {
        source(paste(path, "functions_estimate_rcsdid.R", sep="/"))
        df = df
        treatment = "treatment"
        unit="group"
        time="time"
        estim=estim
        
        df$y_bb = df$y + rnorm(nrow(df))
        

        outcome = "y_bb"
        tres = rcsdid(df=df, treatment=treatment, unit=unit, time=time, outcome=outcome, adj_formula=NULL, estim=estim)$att

        
        return(data.frame(tres, row.names = NULL))
      }
      
      # Save after each case 
        # Prepare output
      if (estim=="all") {
        prepout = data.frame(
          did_bias = mean(res_sim$did)-tre,
          did_sd = sd(res_sim$did-tre),
          did_rmse = sqrt(mean((res_sim$did-tre)^2)),
          rcsdid_bias = mean(res_sim$rcsdid)-tre,
          rcsdid_sd = sd(res_sim$rcsdid-tre),
          rcsdid_rmse = sqrt(mean((res_sim$rcsdid-tre)^2)),
          sdid_bias = mean(res_sim$sdid)-tre,
          sdid_sd = sd(res_sim$sdid-tre),
          sdid_rmse = sqrt(mean((res_sim$sdid-tre)^2))
        )
      }
      if (estim=="rcsdid") {
        prepout = data.frame(
          rcsdid_bias = mean(res_sim$rcsdid)-tre,
          rcsdid_sd = sd(res_sim$rcsdid-tre),
          rcsdid_rmse = sqrt(mean((res_sim$rcsdid-tre)^2))
        )
      }
      stopCluster(cl) # Stop clusters
      
      # Save results
      case = paste(nco, nT, paste(rc_dim ,collapse="_"), "cor_scale", cor_scale, "w", w, "IFE", nbf, sep="_")
      outres[[case]] = prepout
      saveRDS(outres,file=outfile)
    } # end nbf loop
  } # end sample size loop
} # end function 