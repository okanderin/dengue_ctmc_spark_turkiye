## =========================================================
## data_build_once_run.R
## Master script: veri hazırlama + duyarlılık analizi
##
## Çalıştırma: RStudio'da aç, tamamını çalıştır (Ctrl+Shift+Enter)
## Tahmini süre: ~45-60 saat (i7-6700HQ, 16GB RAM)
## =========================================================

cat("\n", strrep("=", 60), "\n")
cat("DATA BUILD PIPELINE — STARTED\n")
cat("Time:", as.character(Sys.time()), "\n")
cat(strrep("=", 60), "\n\n")

## =========================================================
## AŞAMA 1: SSP-bağımsız paylaşımlı veriler (BİR KEZ)
## =========================================================
cat(">>> AŞAMA 1: Paylaşımlı veriler <<<\n\n")

source("R/02_data/build_bias_correction.R")
source("R/02_data/build_population_dataset.R")
source("R/02_data/build_travel_weights.R")
source("R/02_data/build_seasonality_monthly.R")
source("R/02_data/build_parameter_trait.R")
source("R/02_data/build_trait_params_species_files.R")
source("R/02_data/build_sentinel_species.R")

cat("\n>>> Aşama 1 tamamlandı.\n\n")

## =========================================================
## AŞAMA 2: SSP-bağımlı veri hazırlama (HER SSP İÇİN)
## =========================================================
cat(">>> AŞAMA 2: SSP-bağımlı veri hazırlama <<<\n\n")

SSP_LIST <- c("ssp126", "ssp245", "ssp585")

for (ssp in SSP_LIST) {
  cat("\n", strrep("-", 50), "\n")
  cat("  SSP:", ssp, "\n")
  cat(strrep("-", 50), "\n")

  Sys.setenv(SSP_SCENARIO = ssp)

  source("R/02_data/import_raw.R")
  source("R/02_data/clean_interim.R")
  source("R/02_data/build_climate_ssp.R")
  source("R/02_data/build_importation_pressure_monthly.R")
  source("R/02_data/build_rainfall_CSI.R")
}

cat("\n>>> Aşama 2 tamamlandı.\n\n")


## =========================================================
## AŞAMA 3: Duyarlılık analizi (k + η) (HER SSP İÇİN)
##
## Base pipeline'ın tamamlanmış çıktılarına bağlıdır.
## Bu yüzden Aşama 2'den SONRA çalışır.
## =========================================================
cat(">>> AŞAMA 3: Duyarlılık analizi <<<\n\n")

for (ssp in SSP_LIST) {
  cat("\n  Sensitivity:", ssp, "\n")
  Sys.setenv(SSP_SCENARIO = ssp)
  source("R/02_data/run_sensitivity_importation.R")
}

cat("\n>>> Aşama 3 tamamlandı.\n\n")

## =========================================================
## AŞAMA 4: Model çalıştırmagetw (HER SSP İÇİN)
##
## ctmc_spark_monte_carlo.R defines run_ctmc_spark() which:
##   - reads SSP-dependent climate + importation from DIR_PROCESSED_SSP
##   - writes model outputs to DIR_OUTPUT_SSP/simulation
##   - runs n_rep=1000 MC replications with stochastic EIP
##
## Tahmini süre: ~3-5 saat / SSP (i7-6700HQ, 16GB RAM)
## =========================================================
cat(">>> AŞAMA 4: Model çalıştırma <<<\n\n")

for (ssp in SSP_LIST) {
  cat("\n  Model run:", ssp, "\n")
  Sys.setenv(SSP_SCENARIO = ssp)
  source("R/03_models/ctmc_spark_monte_carlo.R")
}

cat("\n>>> Aşama 4 tamamlandı.\n\n")

## =========================================================
## AŞAMA 5: Model duyarlılık analizi (HER SSP İÇİN)
##
## sensitivity_ctmc_mc.R runs 9 scenarios × n_rep=1000 MC:
##   m ∈ {0.50, 0.80, 1.00, 1.20, 2.00}
##   beta_vh ± 20%, ip_days ± 20%
##
## Tahmini süre: ~20-30 saat / SSP (i7-6700HQ, 16GB RAM)
## İlk test: n_rep=50 ile çalıştırarak doğrulayın.
## =========================================================
cat(">>> AŞAMA 5: Model duyarlılık analizi <<<\n\n")

for (ssp in SSP_LIST) {
  cat("\n  Model sensitivity:", ssp, "\n")
  Sys.setenv(SSP_SCENARIO = ssp)
  source("R/03_models/sensitivity_ctmc_mc.R")
}
cat("\n>>> Aşama 5 tamamlandı.\n\n")

cat("\n", strrep("=", 60), "\n")
cat("FULL PIPELINE — FINISHED (Stages 1-5)\n")
cat("Time:", as.character(Sys.time()), "\n")
cat(strrep("=", 60), "\n")
