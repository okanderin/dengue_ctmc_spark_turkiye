for (s in c("ssp126", "ssp245", "ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  
  # Analitik çıktı var mı kontrol et
  analytic_dir <- file.path(here::here("outputs", s, "model_results"))
  has_analytic <- length(list.files(analytic_dir, "ctmc_spark_monthly.*\\.rds$")) > 0
  
  if (has_analytic) {
    source("R/04_results/mc_validation_vs_analytic.R")
  } else {
    message("⚠ ", s, ": analitik çıktı yok — atlanıyor.")
  }
}