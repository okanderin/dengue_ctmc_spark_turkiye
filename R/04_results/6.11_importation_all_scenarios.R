for (s in c("ssp126","ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source("R/02_data/build_importation_pressure_monthly.R")   # tw_path zaten zon
}
source("R/04_results/run_pipeline_3_4.R")   # CTMC'yi 3 SSP için yeniden koşar