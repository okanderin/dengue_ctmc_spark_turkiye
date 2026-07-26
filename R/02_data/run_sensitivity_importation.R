## =========================================================
## run_sensitivity_importation.R (V2 — GBD-weighted)
##
## Sensitivity analysis for importation parameters:
##   Part A: k (climate elasticity) — re-run climate + importation pipeline
##   Part B: η (viremic fraction) — post-processing linear scaling
##
## NOTE: ρ sensitivity is NO LONGER NEEDED because GBD incidence
##   already includes under-reporting correction. GBD uncertainty
##   is captured by low/main/high scenarios in the base pipeline.
##
## References:
##   k: Cheng et al. (2023) eBioMedicine — meta-analysis RR=1.13/°C
##   η: Wilder-Smith et al. (2017) Eurosurveillance; (2019) Lancet
##      Sánchez-Carbonel et al. (2024) J Travel Med
## =========================================================

source("R/01_setup/init.R")

SSP_SCENARIO <- Sys.getenv("SSP_SCENARIO", unset = "ssp245")
cat("\n=============================================\n")
cat("run_sensitivity_importation.R (GBD-weighted, V2)\n")
cat("SSP:", SSP_SCENARIO, "\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("=============================================\n\n")

## ---- Sensitivity configuration ----
K_BASE <- 0.13
K_GRID <- c(0.10, 0.13, 0.20)

ETA_BASE <- 0.25
ETA_GRID <- tibble::tibble(
  label = c("eta_low", "eta_base", "eta_high"),
  eta   = c(0.15, 0.25, 0.33),
  scale = c(0.15, 0.25, 0.33) / ETA_BASE  # linear scaling factor
)

SENS_DIR <- file.path(DIR_OUTPUT_SSP, "sensitivity", "importation")
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## PART A: k sensitivity (re-run climate + importation pipeline)
## =========================================================
cat("\n>>> PART A: k sensitivity <<<\n")
cat("k values:", paste(K_GRID, collapse = ", "), "\n\n")

for (k_val in K_GRID) {
  k_label <- sprintf("k%03d", round(k_val * 100))
  k_dir <- file.path(SENS_DIR, k_label)
  dir.create(k_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("  ", k_label, ": k =", k_val, "\n")
  
  Sys.setenv(K_ELASTICITY = as.character(k_val))
  source("R/02_data/build_climate_ssp.R", local = new.env(parent = globalenv()))
  source("R/02_data/build_importation_pressure_monthly.R",
         local = new.env(parent = globalenv()))
  
  # Copy outputs to sensitivity directory
  for (ext in c(".rds")) {
    for (prefix in c("importation_pressure_monthly_2025_2075",
                     "importation_pressure_yearly_2025_2075",
                     "import_climate_multiplier_monthly")) {
      src <- file.path(DIR_PROCESSED_SSP, paste0(prefix, ext))
      if (file.exists(src)) {
        file.copy(src, file.path(k_dir, paste0(prefix, ext)),
                  overwrite = TRUE)
      }
    }
  }
  
  cat("    Outputs saved to:", k_dir, "\n")
}

# Restore base k
Sys.setenv(K_ELASTICITY = as.character(K_BASE))
cat("  Restoring base k =", K_BASE, "\n")
source("R/02_data/build_climate_ssp.R", local = new.env(parent = globalenv()))
source("R/02_data/build_importation_pressure_monthly.R",
       local = new.env(parent = globalenv()))
cat("  Base pipeline restored.\n\n")

## =========================================================
## PART B: η sensitivity (post-processing)
##
## With GBD-weighted approach, ρ is removed (GBD already corrected).
## η (viremic fraction) is the remaining importation uncertainty.
##
## Since λ_import ∝ η (linearly through the pipeline),
## scaling by η_new/η_base is exact.
##
## η values:
##   0.15 = short viremic window / long travel time
##   0.25 = base (Wilder-Smith et al. 2019)
##   0.33 = long viremic window / short travel time
## =========================================================
cat("\n>>> PART B: η sensitivity <<<\n")
cat("η values:", paste(ETA_GRID$eta, collapse = ", "), "\n\n")

# Load base importation pressure
base_imp_path <- file.path(DIR_PROCESSED_SSP,
                           "importation_pressure_monthly_2025_2075.rds")
if (!file.exists(base_imp_path)) {
  stop("Base importation file not found: ", base_imp_path, call. = FALSE)
}

imp_base <- readRDS(base_imp_path)

for (j in seq_len(nrow(ETA_GRID))) {
  e_label <- ETA_GRID$label[j]
  e_scale <- ETA_GRID$scale[j]
  e_eta   <- ETA_GRID$eta[j]
  e_dir   <- file.path(SENS_DIR, e_label)
  dir.create(e_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("  ", e_label, ": η=", e_eta,
      " scale=", round(e_scale, 3), "\n")
  
  # Scale all lambda columns
  imp_scaled <- imp_base %>%
    dplyr::mutate(
      eta_sens = e_eta,
      scale_factor = e_scale,
      
      lambda_import_per_day = lambda_import_per_day * e_scale,
      lambda_import_window  = lambda_import_window  * e_scale,
      lambda_import         = lambda_import         * e_scale,
      expected_imported_cases_per_month =
        expected_imported_cases_per_month * e_scale,
      
      # Recompute q from scaled lambda
      q_import_month = 1 - exp(-lambda_import_window)
    )
  
  saveRDS(imp_scaled,
          file.path(e_dir, "importation_pressure_monthly_2025_2075.rds"))
  readr::write_csv(imp_scaled,
                   file.path(e_dir, "importation_pressure_monthly_2025_2075.csv"))
  
  # Yearly summary
  imp_y <- imp_scaled %>%
    dplyr::group_by(district_id, province_name, district_name,
                    pi_scenario, year) %>%
    dplyr::summarise(
      lambda_import_year = sum(lambda_import, na.rm = TRUE),
      q_import_year = 1 - prod(1 - q_import_month, na.rm = TRUE),
      M_climate_mean = mean(M_climate, na.rm = TRUE),
      .groups = "drop"
    )
  
  saveRDS(imp_y,
          file.path(e_dir, "importation_pressure_yearly_2025_2075.rds"))
  
  cat("    Saved to:", e_dir, "\n")
}

## =========================================================
## Summary table
## =========================================================
cat("\n>>> Generating sensitivity summary <<<\n")

summary_rows <- list()

# k results
for (k_val in K_GRID) {
  k_label <- sprintf("k%03d", round(k_val * 100))
  k_path <- file.path(SENS_DIR, k_label,
                      "importation_pressure_yearly_2025_2075.rds")
  if (file.exists(k_path)) {
    k_data <- readRDS(k_path)
    for (yr in c(2025, 2050, 2075)) {
      row <- k_data %>%
        dplyr::filter(year == yr, pi_scenario == "main") %>%
        dplyr::summarise(
          lambda_total = sum(lambda_import_year),
          M_climate_mean = mean(M_climate_mean, na.rm = TRUE)
        )
      summary_rows[[length(summary_rows) + 1]] <- tibble::tibble(
        parameter = "k", level = k_label, value = k_val,
        year = yr,
        lambda_total = row$lambda_total,
        M_climate = row$M_climate_mean
      )
    }
  }
}

# η results
for (j in seq_len(nrow(ETA_GRID))) {
  e_label <- ETA_GRID$label[j]
  e_path <- file.path(SENS_DIR, e_label,
                      "importation_pressure_yearly_2025_2075.rds")
  if (file.exists(e_path)) {
    e_data <- readRDS(e_path)
    for (yr in c(2025, 2050, 2075)) {
      row <- e_data %>%
        dplyr::filter(year == yr, pi_scenario == "main") %>%
        dplyr::summarise(
          lambda_total = sum(lambda_import_year),
          M_climate_mean = mean(M_climate_mean, na.rm = TRUE)
        )
      summary_rows[[length(summary_rows) + 1]] <- tibble::tibble(
        parameter = "eta", level = e_label, value = ETA_GRID$eta[j],
        year = yr,
        lambda_total = row$lambda_total,
        M_climate = row$M_climate_mean
      )
    }
  }
}

if (length(summary_rows) > 0) {
  summary_df <- dplyr::bind_rows(summary_rows)
  summary_path <- file.path(SENS_DIR,
                            paste0("sensitivity_summary_", SSP_SCENARIO, ".csv"))
  readr::write_csv(summary_df, summary_path)
  cat("\nSensitivity summary saved:", summary_path, "\n")
  print(summary_df, n = Inf)
}

cat("\n=============================================\n")
cat("Sensitivity analysis COMPLETE\n")
cat("SSP:", SSP_SCENARIO, "\n")
cat("Outputs:", SENS_DIR, "\n")
cat("=============================================\n")