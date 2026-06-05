## =========================================================
## R/02_data/build_rainfall_CSI.R
## Climate Suitability Index (CSI) with rainfall component
##
## SSP-AWARE: Reads from DIR_PROCESSED_SSP, writes to DIR_OUTPUT_SSP
##
## Inputs (SSP-dependent):
##   DIR_PROCESSED_SSP/climate_sentinel_monthly_long_2015_2100.rds
##
## Inputs (shared):
##   DIR_PROCESSED/trait_params_albopictus.csv
##
## Outputs (SSP-dependent):
##   DIR_OUTPUT_SSP/rainfall_CSI/CSI_with_rainfall.csv/.rds
##   DIR_OUTPUT_SSP/rainfall_CSI/CSI_comparison_3vs4.csv
##   DIR_OUTPUT_SSP/rainfall_CSI/m_seasonal_modulation.csv
##   DIR_OUTPUT_SSP/rainfall_CSI/precipitation_monthly_summary.csv
## =========================================================

## ----------------------------------------------------------
## 0) Setup
## ----------------------------------------------------------
if (!exists("DIR_PROCESSED_SSP")) {
  if (file.exists("R/01_setup/init.R")) {
    source("R/01_setup/init.R")
  } else {
    stop("DIR_PROCESSED_SSP not defined. Source init.R first.", call. = FALSE)
  }
}

if (!exists("briere")) {
  if (file.exists("R/03_models/parameter_functions.R")) {
    source("R/03_models/parameter_functions.R")
  } else {
    stop("parameter_functions.R not found.", call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(tibble); library(purrr); library(lubridate)
})

## ----------------------------------------------------------
## 1) Rainfall suitability function: R_suit(P)
## ----------------------------------------------------------
rainfall_suitability <- function(P_mm,
                                  P_min = 20, P_opt_lo = 80,
                                  P_opt_hi = 150, P_max = 400,
                                  R_flush = 0.3) {
  P_mm <- as.numeric(P_mm)
  R <- rep(0, length(P_mm))
  R[P_mm >= P_min & P_mm < P_opt_lo] <-
    (P_mm[P_mm >= P_min & P_mm < P_opt_lo] - P_min) / (P_opt_lo - P_min)
  R[P_mm >= P_opt_lo & P_mm <= P_opt_hi] <- 1.0
  idx_dec <- P_mm > P_opt_hi & P_mm <= P_max
  R[idx_dec] <- 1.0 - (1.0 - R_flush) * (P_mm[idx_dec] - P_opt_hi) / (P_max - P_opt_hi)
  R[P_mm > P_max] <- R_flush
  pmax(pmin(R, 1), 0)
}

## ----------------------------------------------------------
## 2) Extract precipitation
## ----------------------------------------------------------
load_precipitation <- function(year_min = 2025, year_max = 2075) {
  # SSP-specific climate file
  fp_climate <- file.path(DIR_PROCESSED_SSP,
                          "climate_sentinel_monthly_long_2015_2100.rds")
  if (!file.exists(fp_climate)) stop("Missing: ", fp_climate, call. = FALSE)

  df <- readRDS(fp_climate)
  message("  Loaded: ", fp_climate)

  if (!"pr" %in% unique(df$variable)) {
    stop("'pr' not found. Available: ",
         paste(unique(df$variable), collapse = ", "), call. = FALSE)
  }

  df %>%
    filter(variable == "pr", year >= year_min, year <= year_max) %>%
    mutate(
      pr_mm_day = as.numeric(value),
      days_in_month = as.integer(lubridate::days_in_month(
        lubridate::make_date(year, month, 1L))),
      pr_mm = pr_mm_day * days_in_month
    ) %>%
    select(district_id, year, month, pr_mm_day, pr_mm)
}

## ----------------------------------------------------------
## 3) Build extended CSI with rainfall
## ----------------------------------------------------------
build_CSI_with_rainfall <- function(
    year_min = 2025, year_max = 2075,
    out_dir  = file.path(DIR_OUTPUT_SSP, "rainfall_CSI")
) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\n>>> Building CSI with rainfall (", SSP_SCENARIO, ")\n")
  cat("  Output dir:", out_dir, "\n")

  fp_climate    <- file.path(DIR_PROCESSED_SSP,
                             "climate_sentinel_monthly_long_2015_2100.rds")
  fp_traits_alb <- file.path(DIR_PROCESSED, "trait_params_albopictus.csv")

  if (!file.exists(fp_climate))    stop("Missing: ", fp_climate, call. = FALSE)
  if (!file.exists(fp_traits_alb)) stop("Missing: ", fp_traits_alb, call. = FALSE)

  tab_alb <- read_trait_params(fp_traits_alb)

  df_clim <- readRDS(fp_climate)
  if (all(c("variable", "value") %in% names(df_clim))) {
    df_clim <- standardize_climate_cols(df_clim,
                                        key = c("district_id", "year", "month"))
  }
  df_clim <- df_clim %>% filter(year >= year_min, year <= year_max)

  df_precip <- load_precipitation(year_min, year_max)

  df <- df_clim %>%
    left_join(df_precip, by = c("district_id", "year", "month"))

  if (any(is.na(df$pr_mm))) {
    df$pr_mm[is.na(df$pr_mm)] <- 0
    df$pr_mm_day[is.na(df$pr_mm_day)] <- 0
  }

  r_a   <- tab_alb %>% filter(trait == "a_biting")
  r_eip <- tab_alb %>% filter(trait == "eip_dev_rate")
  r_lf  <- tab_alb %>% filter(trait == "lifespan_temp")
  lf_floor <- 0.25

  df_csi <- df %>%
    mutate(
      a_raw  = briere(temp_c, r_a$c, r_a$T0, r_a$Tm),
      lf_raw = pmax(quadratic_unimodal(temp_c, r_lf$c, r_lf$T0, r_lf$Tm),
                    lf_floor),
      eip_dr = briere(temp_c, r_eip$c, r_eip$T0, r_eip$Tm),
      R_suit = rainfall_suitability(pr_mm)
    ) %>%
    mutate(
      a_norm   = a_raw  / max(a_raw,  na.rm = TRUE),
      lf_norm  = lf_raw / max(lf_raw, na.rm = TRUE),
      eip_norm = eip_dr / max(eip_dr, na.rm = TRUE),
      R_norm   = R_suit,
      CSI_3 = (a_norm + lf_norm + eip_norm) / 3,
      CSI_4 = (a_norm + lf_norm + eip_norm + R_norm) / 4
    )

  write_csv(df_csi, file.path(out_dir, "CSI_with_rainfall.csv"))
  saveRDS(df_csi, file.path(out_dir, "CSI_with_rainfall.rds"))

  csi_comparison <- df_csi %>%
    group_by(district_id, month) %>%
    summarise(CSI_3_mean = mean(CSI_3), CSI_4_mean = mean(CSI_4),
              R_suit_mean = mean(R_suit), pr_mm_mean = mean(pr_mm),
              delta_CSI = mean(CSI_4 - CSI_3), .groups = "drop")
  write_csv(csi_comparison, file.path(out_dir, "CSI_comparison_3vs4.csv"))

  m_mod <- df_csi %>%
    group_by(district_id, month) %>%
    summarise(R_suit_mean = mean(R_suit), pr_mm_mean = mean(pr_mm),
              .groups = "drop") %>%
    mutate(m_eff_factor = R_suit_mean)
  write_csv(m_mod, file.path(out_dir, "m_seasonal_modulation.csv"))

  precip_sum <- df_csi %>%
    group_by(district_id, month) %>%
    summarise(pr_mm_mean = mean(pr_mm), pr_mm_sd = sd(pr_mm),
              .groups = "drop")
  write_csv(precip_sum, file.path(out_dir, "precipitation_monthly_summary.csv"))

  cat("  CSI outputs saved to:", out_dir, "\n")
  invisible(df_csi)
}

## ----------------------------------------------------------
## Entrypoint
## ----------------------------------------------------------
if (identical(environment(), globalenv())) {
  rainfall_csi_results <- build_CSI_with_rainfall()
}
message("build_rainfall_CSI.R loaded (", SSP_SCENARIO, ")")
