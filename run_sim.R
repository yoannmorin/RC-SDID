###########################################
###########################################
###########################################
#####                                 #####
#####                                 #####
##### Output simulation table results #####
#####                                 #####
#####                                 #####
###########################################
###########################################
###########################################

 
###### Start benchmark
start_time <- Sys.time()


##### Set seed
set.seed(123)


##### Path to directory
path="yourpathhere"

setwd(path)


##### Packages and functions
source("functions_gen_rng_var.R")
source("functions_simul_rcsdid.R")


##### Set parameters common to all simulations
nrep = 1000
estim="all"
tre = 0.3


##### Create random variables to be held constant across all simulations
functions_gen_rng_var(path)







###########################
###########################
#####                 #####
##### Run simulations #####
#####                 #####
###########################
###########################

##### Run table scale parameter ##### 
vnbf = c(1)
outfile = "tab_scale.rds"
samp_size = list(list("30", "30", "scalenor", "1", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "2", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "4", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "6", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "8", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "10", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "15", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "20", "100", "02", "0.2")
)

functions_simul_rcsdid(nrep, vnbf, samp_size, path, tre, estim, outfile)






##### Run table number of IFE #####
samp_size = list(list("30", "30", "scalenor", "10", "100", "02", "0.2"))
outfile = "tab_ife.rds"
vnbf = c(0,1,2,3,4)

functions_simul_rcsdid(nrep, vnbf, samp_size, path, tre, estim, outfile)






##### Run table correlation individual FE and scale parameter ##### 
vnbf = c(1)
outfile = "tab_cors.rds"
samp_size = list(list("30", "30", "scalenor", "10", "100", "0", "0.2"),
                 list("30", "30", "scalenor", "10", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "10", "100", "05", "0.2"),
                 list("30", "30", "scalenor", "10", "100", "08", "0.2"),
                 list("30", "30", "scalenor", "10", "100", "1", "0.2")
)

functions_simul_rcsdid(nrep, vnbf, samp_size, path, tre, estim, outfile)





##### Run table correlation treatment FE/IFE ##### 
vnbf = c(1)
outfile = "tab_w.rds"
samp_size = list(list("30", "30", "scalenor", "10", "100", "02", "1"),
                 list("30", "30", "scalenor", "10", "100", "02", "0.8"),
                 list("30", "30", "scalenor", "10", "100", "02", "0.6"),
                 list("30", "30", "scalenor", "10", "100", "02", "0.4"),
                 list("30", "30", "scalenor", "10", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "10", "100", "02", "0")
                 )

functions_simul_rcsdid(nrep, vnbf, samp_size, path, tre, estim, outfile)




##### Run table size ##### 
vnbf = c(1)
outfile = "tab_size.rds"
samp_size = list(list("30", "30", "scalenor", "10", "100", "02", "0.2"),
                 list("15", "30", "scalenor", "10", "100", "02", "0.2"),
                 list("30", "15", "scalenor", "10", "100", "02", "0.2"),
                 list("15", "15", "scalenor", "10", "100", "02", "0.2"),
                 list("30", "30", "scalenor", "10", "50", "02", "0.2"),
                 list("15", "30", "scalenor", "10", "50", "02", "0.2"),
                 list("30", "15", "scalenor", "10", "50", "02", "0.2"),
                 list("15", "15", "scalenor", "10", "50", "02", "0.2")
)

functions_simul_rcsdid(nrep, vnbf, samp_size, path, tre, estim, outfile)







###### End benchmark
end_time <- Sys.time()
end_time - start_time
#  4.573785 hours using 16 threads on a Ryzen 7 7700X.







###########################
###########################
#####                 #####
#####  Output tables  #####
#####                 #####
###########################
###########################

##### Packages
library("dplyr")

##### Working directory
setwd(path)



##### Import results
tab_scale = readRDS("tab_scale.rds")
tab_ife = readRDS("tab_ife.rds")
tab_cors = readRDS("tab_cors.rds")
tab_w = readRDS("tab_w.rds")
tab_size = readRDS("tab_size.rds")


##### Final tables
fintab = function(df){
  fdf = NULL
  for (j in 1:length(df)) {
    namerow = names(df)[[j]]
    adf = data.frame(df[[j]])
    rownames(adf) = namerow
    fdf = rbind(fdf, adf)
  }
  rdf = fdf[,c(1,4,7,2,5,8,3,6,9)]
  return(rdf)
}

ftab_scale = fintab(tab_scale)
ftab_ife = fintab(tab_ife)
ftab_cors = fintab(tab_cors)
ftab_w = fintab(tab_w)
ftab_size = fintab(tab_size)




# Output tex
ftab_scale %>%
  knitr::kable(format = 'latex', booktabs = TRUE) %>%
  writeLines('sim_tab_scale.tex')

ftab_ife %>%
  knitr::kable(format = 'latex', booktabs = TRUE) %>%
  writeLines('sim_tab_ife.tex')

ftab_cors %>%
  knitr::kable(format = 'latex', booktabs = TRUE) %>%
  writeLines('sim_tab_cors.tex')

ftab_w %>%
  knitr::kable(format = 'latex', booktabs = TRUE) %>%
  writeLines('sim_tab_w.tex')

ftab_size %>%
  knitr::kable(format = 'latex', booktabs = TRUE) %>%
  writeLines('sim_tab_size.tex')
