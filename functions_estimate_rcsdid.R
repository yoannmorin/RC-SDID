## -----------------------------------------------------------------------------
## functions_estimate_rcsdid.R
## Yoann Morin
## 
## This creates R functions to compute RC-SDID estimate.
## Options also allow to compute other estimates for comparison.
## -----------------------------------------------------------------------------

# Packages
library("fixest")
library("dplyr")
library("synthdid")

#' Compute RC-SDID estimate
#' 
#' @param df - Dataframe 
#' @param treatment - binary treatment variable (1 if i is treated at time t, 0 else)
#' @param unit - ID variable for the unit level at the level the data will be aggregated at
#' @param time - time variable (also used to aggregate data)
#' @param outcome -  name or the outcome variable if there are no covariates, "y_adj" if there are covariates
#' @param adj_formula - fixest formula to adjust for covariates (individual and time fixed effects can be at a different level than the one used to aggregate data)
#' @param estim - "rcsdid" to compute rcsdid estimate (the default), "all" to also compute DID and sdid without RC adjustment.
rcsdid <- function(df, treatment, unit, time, outcome, adj_formula=NULL, estim=NULL){
  
  if(!is.null(adj_formula)) { 
  print("Adjusting for covariates")
  
  # Compute formualas for ATT
  ffe = gsub(".*\\|","",adj_formula)[3]
  fm = as.formula(paste(paste(outcome, treatment, sep="~"), ffe, sep="|"))
  fmd = as.formula(paste(paste(gsub("\\|.*","",adj_formula)[2], paste(gsub("\\|.*","",adj_formula)[3], treatment, sep="+"), sep="~"), ffe, sep="|"))
  
  # Adjust for covariates
  tdf = df[df[treatment] ==0, ] 
  reg_adj =  feols(adj_formula, data = tdf)
  
  tform = paste("~", gsub("\\|.*","",gsub(".*~","",adj_formula)[3]), sep="")
  xform =  as.formula(tform)
  xmat = model.matrix(xform, df)[,-1, drop=FALSE]
  xvars = colnames(xmat)
  
  beta = coef(reg_adj)[seq_along(xvars)]
  xeff = as.numeric(xmat %*% beta)
  
  var_dep = gsub(".*~","",adj_formula)[2]
  
  df$y_adj = df[[var_dep]] - xeff 
  }
  else {
  print("No covariates provided.")
  # Compute formualas for ATT
  fm = as.formula(paste(outcome, paste(treatment, paste(unit, time, sep="+"), sep="|"), sep="~"))
  fmd = fm
  }  
  
  ### Compute weights
  print("Computing weights.")
  # Aggregate data
  agg_vars = c(unit, time)
  ovars =  c(unit, time, treatment, outcome)
  
  adf = df %>%
    group_by(across(all_of(agg_vars))) %>%
    select(all_of(ovars)) %>%
    summarise_all(list(mean))
  
  # Compute weight using the SDID method in the synthdid package
  adf = as.data.frame(adf)
  adf[[treatment]] = as.integer(adf[[treatment]])
  msetup = panel.matrices(adf, unit=unit, time=time, outcome=outcome, treatment=treatment)
  
  sdid = synthdid_estimate(msetup$Y, msetup$N0, msetup$T0)
  
  # Weights
  womega = data.frame(unit_weights = head(rownames(attr(sdid, 'setup')$Y),attr(sdid, 'setup')$N0) , omega_wei = attr(sdid, 'weights')$omega)
  wlambda = data.frame(time_weights = head(colnames(attr(sdid, 'setup')$Y),attr(sdid, 'setup')$T0), lambda_wei = wlambda <- attr(sdid, 'weights')$lambda)
  
  
  # All weights in database
  utreat = adf %>%
    select(all_of(c(unit, treatment))) %>%
    filter_all(any_vars(c(treatment)==1)) %>%
    distinct() %>%
    select(all_of(c(unit)))
  utreat = as.character(utreat)
  
  fullomega = rbind(womega, data.frame(unit_weights =utreat, omega_wei = 1 / (nrow(attr(sdid, 'setup')$Y) - attr(sdid, 'setup')$N0)))
  fulllambda = rbind(wlambda, data.frame(time_weights = tail(colnames(attr(sdid, 'setup')$Y),ncol(attr(sdid, 'setup')$Y) - attr(sdid, 'setup')$T0), lambda_wei = 1 / (ncol(attr(sdid, 'setup')$Y) -  attr(sdid, 'setup')$T0)))
  
  df = merge(df, fullomega, by.x=c(unit), by.y=c("unit_weights"))
  df = merge(df, fulllambda, by=c(time), by.y=c("time_weights"))
  
  df$weight_sdid = df$omega_wei*df$lambda_wei
  
  
  
  ### Compute weights specific for RC SDID
  df$treat_fe = ifelse(df[unit]==utreat, 1, 0)
  
  ## For control group
  # Number of observations in group k in period t
  wco = df %>%
    filter(treat_fe == 0) %>%
    group_by(across(all_of(agg_vars))) %>%
    summarise(nbobs_kt = n())
  
  # Compute rc weights
  wco$weight_rc = (1 / wco$nbobs_kt)
  
  ## For treated group
  # Number of treated groups
  wtr = df %>%
    filter(treat_fe == 1) %>%
    group_by(across(all_of(agg_vars))) %>%
    summarise(nbobs_kt = n())
  
  # Compute rc weights
  wtr$weight_rc = (1 / wtr$nbobs_kt)
  
  ## Merge results
  wgroup = rbind(wco, wtr)
  df = merge(df, wgroup, by=c(unit, time))
  
  ## Compute rcsdid weights
  df$weight_rcsdid = df$weight_sdid*df$weight_rc
  
  

  ### Compute ATT
  print("Computing ATT.")
 
  
  # rcsdid 
  func_res_rcsdid = feols(df, fm, weights = df$weight_rcsdid)$coefficients[treatment]
  
  
  if (estim=="all") {
  # sdid without RC adjustment
  func_res_sdid = feols(df, fm, weights = df$weight_sdid)$coefficients[treatment]
  
  # did
  func_res_did = feols(df, fmd)$coefficients[treatment]
  }
  
  # Results in a table
  estimates = data.frame(do.call(cbind, mget(ls(pattern = 'func_res_'))))
  colnames(estimates) = substr(colnames(estimates), 10, 30)
  
  
  
  ### Output
  out_function = c()
  out_function$att = estimates
  out_function$unit_weights = womega
  out_function$time_weights = wlambda
  out_function$rcsdid_weights = wgroup
  
  return(out_function)
}