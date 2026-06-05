for (s in c("ssp126", "ssp245", "ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source("R/04_results/mc_validation_vs_analytic.R")
}