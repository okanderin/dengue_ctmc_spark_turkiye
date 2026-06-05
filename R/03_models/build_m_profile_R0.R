# ==========================================================
# build_m_profile_R0.R
# Profile analysis of m parameter, m_critical calculation,
# and R0 distribution reporting
# SSP-AWARE: reads climate from DIR_PROCESSED_SSP, writes to DIR_OUTPUT_SSP
#
# PURPOSE:
#   1) For each district-month, compute R0 as a function of m
#   2) Find m_critical: the m value where R0 = 1 (subcritical/supercritical boundary)
#   3) Report full R0 distribution under m ~ Uniform(0.1, 2.0)
#   4) Profile-likelihood style sweep: P_est(m) for each district across m grid
#   5) Bayesian posterior approximation for m given external prevalence constraints
#
# INPUTS:
#   - climate_sentinel_monthly_long_2015_2100.rds  (temperature, RH)
#   - trait_params_aegypti.csv / trait_params_albopictus.csv
#   - sentinel_species.csv
#
# OUTPUTS:
#   - m_critical_by_district_month.csv      (m_critical per district-month-year)
#   - R0_distribution_summary.csv           (R0 stats under m ~ Uniform)
#   - m_profile_P_est.csv                   (P_est vs m sweep per district)
#   - R0_monthly_full_grid.rds              (full m x district x month grid)
#
# RUN:
#   source("R/01_setup/init.R")
#   source("R/03_models/parameter_functions.R")
#   source("R/03_models/build_m_profile_R0.R")
# ==========================================================

## ----------------------------------------------------------
## 0) Setup — load only if not already loaded by caller
## ----------------------------------------------------------
if (!exists("DIR_PROCESSED")) {
  if (file.exists("R/01_setup/init.R")) {
    source("R/01_setup/init.R")
  } else {
    stop("DIR_PROCESSED not defined. Source init.R first or define paths.", call. = FALSE)
  }
}

if (!exists("make_lambda_local")) {
  if (file.exists("R/03_models/parameter_functions.R")) {
    source("R/03_models/parameter_functions.R")
  } else {
    stop("parameter_functions.R not found. It must define make_lambda_local().", call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
})

## ----------------------------------------------------------
## 1) Helper: analytic R0 from climate + parameters
##    R0 = (m * a(T)^2 * beta_vh * beta_hv * exp(-mu_v * EIP)) / (mu_v * gamma)
##
##    Returns R0 for given (T, RH, m, beta_vh, beta_hv, gamma, tab_traits)
## ----------------------------------------------------------
compute_R0_analytic <- function(T, RH, m, beta_vh, beta_hv, gamma_day, tab_traits) {
  aT   <- a_of_T(T, tab_traits)
  muv  <- mu_v_of_TRH(T, RH, tab_traits)
  eip  <- eip_mean_of_T(T, tab_traits)

  # Outside thermal range: no transmission

  if (!is.finite(eip) || eip <= 0 || aT <= 0 || muv <= 0) {
    return(0)
  }

  surv <- exp(-muv * eip)
  R0   <- (m * aT^2 * beta_vh * beta_hv * surv) / (muv * gamma_day)
  max(R0, 0)
}

## ----------------------------------------------------------
## 2) Helper: m_critical — solve R0(m_crit) = 1 analytically
##    From R0 formula: m_crit = (mu_v * gamma) / (a^2 * beta_vh * beta_hv * S)
##    where S = exp(-mu_v * EIP)
## ----------------------------------------------------------
compute_m_critical <- function(T, RH, beta_vh, beta_hv, gamma_day, tab_traits) {
  aT   <- a_of_T(T, tab_traits)
  muv  <- mu_v_of_TRH(T, RH, tab_traits)
  eip  <- eip_mean_of_T(T, tab_traits)

  if (!is.finite(eip) || eip <= 0 || aT <= 0 || muv <= 0) {
    return(Inf)  # no transmission possible at any m
  }

  surv <- exp(-muv * eip)
  denom <- aT^2 * beta_vh * beta_hv * surv

  if (denom <= 0) return(Inf)

  m_crit <- (muv * gamma_day) / denom
  m_crit
}

## ----------------------------------------------------------
## 3) Helper: P_est from birth-death theory
##    For R0 > 1:  P_est = 1 - 1/R0
##    For R0 <= 1: P_est = 0  (deterministic; stochastic adds tiny deviations)
## ----------------------------------------------------------
P_est_analytic <- function(R0) {
  ifelse(R0 > 1, 1 - 1/R0, 0)
}

## ----------------------------------------------------------
## 4) Main function: build m profile analysis
## ----------------------------------------------------------
build_m_profile <- function(
    year_min  = 2025,
    year_max  = 2075,
    beta_vh   = 0.30,
    beta_hv   = 0.33,
    infectious_period_days = 5,
    m_grid    = c(seq(0.05, 0.50, by = 0.05),
                  seq(0.60, 1.00, by = 0.10),
                  seq(1.25, 3.00, by = 0.25)),
    m_prior_range = c(0.1, 2.0),   # Uniform prior for R0 distribution
    n_m_prior = 500L,               # samples from m prior for distribution
    out_dir   = file.path(DIR_OUTPUT_SSP, "m_profile")
) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gamma_day <- 1 / infectious_period_days

  # ---- Load inputs (SSP-dependent climate, SSP-independent traits) ----
  fp_climate <- file.path(DIR_PROCESSED_SSP, "climate_sentinel_monthly_long_2015_2100.rds")
  fp_traits_aeg <- file.path(DIR_PROCESSED, "trait_params_aegypti.csv")
  fp_traits_alb <- file.path(DIR_PROCESSED, "trait_params_albopictus.csv")
  fp_species    <- file.path(DIR_PROCESSED, "sentinel_species.csv")

  if (!file.exists(fp_climate))    stop("Missing: ", fp_climate, call. = FALSE)
  if (!file.exists(fp_traits_aeg)) stop("Missing: ", fp_traits_aeg, call. = FALSE)
  if (!file.exists(fp_traits_alb)) stop("Missing: ", fp_traits_alb, call. = FALSE)

  tab_aeg <- read_trait_params(fp_traits_aeg)
  tab_alb <- read_trait_params(fp_traits_alb)

  # Load and standardize climate
  df_clim <- readRDS(fp_climate)
  # Standardize: ensure temp_c and rh columns
  if (all(c("variable", "value") %in% names(df_clim))) {
    df_clim <- standardize_climate_cols(df_clim, key = c("district_id", "year", "month"))
  }

  df_clim <- df_clim %>%
    filter(year >= year_min, year <= year_max)

  # Species mapping
  if (file.exists(fp_species)) {
    sp_map <- read_csv(fp_species, show_col_types = FALSE) %>%
      mutate(species = tolower(species))
  } else {
    ids_aeg <- c("TUR.10.4_1", "TUR.81.6_1")
    sp_map <- tibble(
      district_id = unique(df_clim$district_id),
      species = ifelse(district_id %in% ids_aeg, "aegypti", "albopictus")
    )
  }

  df_clim <- df_clim %>% left_join(sp_map, by = "district_id")

  message("=== m Profile Analysis ===")
  message("Districts: ", n_distinct(df_clim$district_id))
  message("Year range: ", year_min, "–", year_max)
  message("m grid points: ", length(m_grid))
  message("m prior: Uniform(", m_prior_range[1], ", ", m_prior_range[2], ")")

  # ================================================================
  # PART A: m_critical for each district-month-year
  # ================================================================
  message("\n[A] Computing m_critical for each district-month-year...")

  df_mcrit <- df_clim %>%
    mutate(
      tab_traits_name = species,
      m_critical = pmap_dbl(
        list(temp_c, rh, species),
        function(T, RH, sp) {
          tab <- if (sp == "aegypti") tab_aeg else tab_alb
          compute_m_critical(T, RH, beta_vh, beta_hv, gamma_day, tab)
        }
      ),
      R0_at_m050 = pmap_dbl(
        list(temp_c, rh, species),
        function(T, RH, sp) {
          tab <- if (sp == "aegypti") tab_aeg else tab_alb
          compute_R0_analytic(T, RH, m = 0.50, beta_vh, beta_hv, gamma_day, tab)
        }
      ),
      R0_at_m100 = pmap_dbl(
        list(temp_c, rh, species),
        function(T, RH, sp) {
          tab <- if (sp == "aegypti") tab_aeg else tab_alb
          compute_R0_analytic(T, RH, m = 1.00, beta_vh, beta_hv, gamma_day, tab)
        }
      ),
      # Is m_critical within plausible range?
      m_crit_plausible = m_critical >= m_prior_range[1] & m_critical <= m_prior_range[2],
      # At base m=0.5, is system supercritical?
      supercritical_m050 = R0_at_m050 > 1,
      supercritical_m100 = R0_at_m100 > 1
    )

  # Save m_critical results
  fp_mcrit <- file.path(out_dir, "m_critical_by_district_month.csv")
  write_csv(df_mcrit, fp_mcrit)
  message("  Saved: ", fp_mcrit)

  # Summary statistics per district
  mcrit_summary <- df_mcrit %>%
    group_by(district_id, species) %>%
    summarise(
      m_crit_min       = min(m_critical, na.rm = TRUE),
      m_crit_median    = median(m_critical[is.finite(m_critical)], na.rm = TRUE),
      m_crit_mean      = mean(m_critical[is.finite(m_critical)], na.rm = TRUE),
      m_crit_p25       = quantile(m_critical[is.finite(m_critical)], 0.25, na.rm = TRUE),
      m_crit_p75       = quantile(m_critical[is.finite(m_critical)], 0.75, na.rm = TRUE),
      n_months_total   = n(),
      n_months_supercrit_m050 = sum(supercritical_m050, na.rm = TRUE),
      n_months_supercrit_m100 = sum(supercritical_m100, na.rm = TRUE),
      n_months_m_crit_plausible = sum(m_crit_plausible, na.rm = TRUE),
      pct_supercrit_m050 = round(mean(supercritical_m050, na.rm = TRUE) * 100, 2),
      pct_supercrit_m100 = round(mean(supercritical_m100, na.rm = TRUE) * 100, 2),
      # Best summer month (lowest m_critical = easiest to cross R0=1)
      best_month = {
        finite_idx <- which(is.finite(m_critical))
        if (length(finite_idx) > 0) month[finite_idx[which.min(m_critical[finite_idx])]] else NA_integer_
      },
      m_crit_best_month = min(m_critical[is.finite(m_critical)], na.rm = TRUE),
      .groups = "drop"
    )

  fp_mcrit_sum <- file.path(out_dir, "m_critical_summary.csv")
  write_csv(mcrit_summary, fp_mcrit_sum)
  message("  Saved: ", fp_mcrit_sum)

  # ================================================================
  # PART B: R0 distribution under m ~ Uniform(0.1, 2.0)
  # ================================================================
  message("\n[B] Computing R0 distribution under m prior...")

  # For each district, compute R0 at peak transmission month (lowest m_critical)
  # across the m prior samples
  set.seed(42)
  m_samples <- runif(n_m_prior, m_prior_range[1], m_prior_range[2])

  # Get representative "peak month" climate for each district
  # (month with lowest m_critical, averaged across years)
  peak_climate <- df_mcrit %>%
    filter(is.finite(m_critical)) %>%
    group_by(district_id, species, month) %>%
    summarise(
      mean_temp = mean(temp_c, na.rm = TRUE),
      mean_rh   = mean(rh, na.rm = TRUE),
      mean_mcrit = mean(m_critical, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(district_id, species) %>%
    slice_min(mean_mcrit, n = 1, with_ties = FALSE) %>%
    ungroup()

  R0_dist <- peak_climate %>%
    rowwise() %>%
    mutate(
      R0_draws = list({
        tab <- if (species == "aegypti") tab_aeg else tab_alb
        sapply(m_samples, function(m_val) {
          compute_R0_analytic(mean_temp, mean_rh, m_val, beta_vh, beta_hv, gamma_day, tab)
        })
      }),
      R0_mean   = mean(unlist(R0_draws)),
      R0_median = median(unlist(R0_draws)),
      R0_p2_5   = quantile(unlist(R0_draws), 0.025),
      R0_p25    = quantile(unlist(R0_draws), 0.25),
      R0_p75    = quantile(unlist(R0_draws), 0.75),
      R0_p97_5  = quantile(unlist(R0_draws), 0.975),
      R0_min    = min(unlist(R0_draws)),
      R0_max    = max(unlist(R0_draws)),
      prob_R0_gt1 = mean(unlist(R0_draws) > 1),
      peak_month_label = month.abb[month]
    ) %>%
    ungroup() %>%
    select(-R0_draws)

  fp_R0dist <- file.path(out_dir, "R0_distribution_summary.csv")
  write_csv(R0_dist, fp_R0dist)
  message("  Saved: ", fp_R0dist)

  # ================================================================
  # PART C: Profile sweep — P_est(m) for each district (peak month)
  # ================================================================
  message("\n[C] Computing P_est profile across m grid...")

  profile_results <- peak_climate %>%
    rowwise() %>%
    mutate(
      profile = list({
        tab <- if (species == "aegypti") tab_aeg else tab_alb
        tibble(
          m_value = m_grid,
          R0 = sapply(m_grid, function(m_val) {
            compute_R0_analytic(mean_temp, mean_rh, m_val, beta_vh, beta_hv, gamma_day, tab)
          }),
          P_est = P_est_analytic(R0)
        )
      })
    ) %>%
    ungroup() %>%
    select(district_id, species, month, mean_temp, mean_rh, mean_mcrit, profile) %>%
    unnest(profile)

  fp_profile <- file.path(out_dir, "m_profile_P_est.csv")
  write_csv(profile_results, fp_profile)
  message("  Saved: ", fp_profile)

  # ================================================================
  # PART D: Full m × district × month grid (for Rmd integration)
  # ================================================================
  message("\n[D] Building full monthly R0 grid (may take a minute)...")

  # Use a coarser m grid for the full dataset to keep memory manageable
  m_grid_coarse <- c(seq(0.1, 0.5, by = 0.1), seq(0.75, 2.0, by = 0.25))

  # Compute for each unique district-month (averaged across years)
  monthly_avg <- df_clim %>%
    group_by(district_id, species, month) %>%
    summarise(
      mean_temp = mean(temp_c, na.rm = TRUE),
      mean_rh   = mean(rh, na.rm = TRUE),
      .groups = "drop"
    )

  R0_grid_full <- monthly_avg %>%
    rowwise() %>%
    mutate(
      R0_by_m = list({
        tab <- if (species == "aegypti") tab_aeg else tab_alb
        tibble(
          m_value = m_grid_coarse,
          R0 = sapply(m_grid_coarse, function(m_val) {
            compute_R0_analytic(mean_temp, mean_rh, m_val, beta_vh, beta_hv, gamma_day, tab)
          }),
          P_est = P_est_analytic(R0),
          supercritical = R0 > 1
        )
      })
    ) %>%
    ungroup() %>%
    unnest(R0_by_m)

  fp_grid <- file.path(out_dir, "R0_monthly_full_grid.rds")
  saveRDS(R0_grid_full, fp_grid)

  fp_grid_csv <- file.path(out_dir, "R0_monthly_full_grid.csv")
  write_csv(R0_grid_full, fp_grid_csv)
  message("  Saved: ", fp_grid, " + .csv")

  # ================================================================
  # PART E: Seasonal m_critical profile (monthly pattern)
  # ================================================================
  message("\n[E] Computing seasonal m_critical pattern...")

  seasonal_mcrit <- df_mcrit %>%
    filter(is.finite(m_critical)) %>%
    group_by(district_id, species, month) %>%
    summarise(
      m_crit_mean   = mean(m_critical, na.rm = TRUE),
      m_crit_median = median(m_critical, na.rm = TRUE),
      m_crit_p10    = quantile(m_critical, 0.10, na.rm = TRUE),
      m_crit_p90    = quantile(m_critical, 0.90, na.rm = TRUE),
      mean_temp     = mean(temp_c, na.rm = TRUE),
      mean_rh       = mean(rh, na.rm = TRUE),
      .groups = "drop"
    )

  fp_seasonal <- file.path(out_dir, "m_critical_seasonal.csv")
  write_csv(seasonal_mcrit, fp_seasonal)
  message("  Saved: ", fp_seasonal)

  # ================================================================
  # Console report
  # ================================================================
  message("\n", strrep("=", 60))
  message("m PROFILE ANALYSIS COMPLETE")
  message(strrep("=", 60))

  message("\n--- m_critical Summary (peak transmission month) ---")
  for (i in seq_len(nrow(mcrit_summary))) {
    row <- mcrit_summary[i, ]
    message(sprintf(
      "  %s (%s): m_crit_min = %.2f | m_crit_median = %.2f | supercrit @ m=0.5: %s | supercrit @ m=1.0: %s",
      row$district_id, row$species,
      row$m_crit_min, row$m_crit_median,
      paste0(row$pct_supercrit_m050, "%"),
      paste0(row$pct_supercrit_m100, "%")
    ))
  }

  message("\n--- R0 Distribution (peak month, m ~ Uniform) ---")
  for (i in seq_len(nrow(R0_dist))) {
    row <- R0_dist[i, ]
    message(sprintf(
      "  %s [%s]: R0 median = %.3f [%.3f–%.3f] | P(R0>1) = %.1f%%",
      row$district_id, row$peak_month_label,
      row$R0_median, row$R0_p2_5, row$R0_p97_5,
      row$prob_R0_gt1 * 100
    ))
  }

  message("\n--- Key Interpretation ---")
  message("  m_critical = m value where R0 crosses 1.0")
  message("  If m_critical < 0.5 : supercritical at base scenario")
  message("  If m_critical ∈ [0.5, 1.0] : sub-to-supercritical transition within OAT range")
  message("  If m_critical > 2.0 : subcritical across entire plausible m range")

  message("\nOutputs: ", out_dir)

  invisible(list(
    m_critical   = df_mcrit,
    mcrit_summary = mcrit_summary,
    R0_dist      = R0_dist,
    profile      = profile_results,
    R0_grid      = R0_grid_full,
    seasonal     = seasonal_mcrit
  ))
}

## ----------------------------------------------------------
## Entrypoint
## ----------------------------------------------------------
if (identical(environment(), globalenv())) {
  m_profile_results <- build_m_profile(
    year_min  = 2025,
    year_max  = 2075,
    beta_vh   = 0.30,
    beta_hv   = 0.33,
    infectious_period_days = 5,
    out_dir   = file.path(DIR_OUTPUT_SSP, "m_profile")
  )
}

message("✓ build_m_profile_R0.R loaded (SSP-aware paths)")
