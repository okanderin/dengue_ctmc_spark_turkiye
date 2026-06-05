# Üç SSP birden:
for (s in c("ssp126", "ssp245", "ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source("R/validation_two_stage.R")
}