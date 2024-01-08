# RC-SDID
Replication Files for Morin Yoann, "Synthetic Difference in Differences for Repeated Cross-Sectional Data"

## Dependencies 
Install the following packages: dplyr, doParallel, foreach, doRNG, synthdid, fixest, rngtools.
install.packages("dplyr")
install.packages("doParallel")
install.packages("foreach")
install.packages("doRNG")
install.packages("fixest")
install.packages("rngtools")
install.packages("devtools")
devtools::install_github("synth-inference/synthdid")

## Replication files for simulations
functions_estimate_rcsdid.R : estimate RCSDID (also estimate SDID and DID)
functions_dgp_rcsdid.R : generate the samples
functions_simul_rcsdid.R : perform the simulations
functions_gen_rng_var.R : draw factors, loadings, covariates and repeated cross sections size only once for all simulations
run_sim.R : simulate then output all tables
