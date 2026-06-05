source("R/01_setup/init.R")

for (s in c("ssp126", "ssp245", "ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source("R/03_models/ctmc_spark.R")
}


source("R/01_setup/init.R")

for (s in c("ssp126", "ssp245", "ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source("R/02_data/build_importation_pressure_monthly.R")
}



