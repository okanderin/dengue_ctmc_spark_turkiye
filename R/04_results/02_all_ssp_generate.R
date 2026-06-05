
for (s in c("ssp126","ssp245","ssp585"))  {
    Sys.setenv(SSP_SCENARIO = s)
    source("R/04_results/01_generate_ssp_outputs.R")
  }



