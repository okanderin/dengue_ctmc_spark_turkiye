## =========================================================
## Project : r_project_tez
## Script  : R/02_data/build_importation_pressure_monthly.R
## Purpose : Monthly importation pressure dataset (sentinel districts)
##           SSP-aware: M_climate varies by scenario.
##
## VERSION 3 — GBD-weighted, S_m relocated to V_in(t)
##
## KEY CHANGES FROM v2:
##
## (1) S_m moved from daily_rate to V_in_day
##     BEFORE: daily_rate = rate_pp_year × S_m × 12 / 365 × M_climate
##             V_in_day   = V_TR × N_total / 365
##     AFTER:  daily_rate = rate_pp_year / 365 × M_climate
##             V_in_day   = V_TR × N_total × S_m / 365
##
##     Rationale: S_m is derived from monthly visitor departure counts
##     (TÜİK/EGM passport records), measuring the seasonal variation
##     in international travel volume to Türkiye. It modulates the
##     number of arriving travelers (exposure population), not the
##     per-capita infection risk in source countries (γ_GBD).
##     The GBD incidence rate is an annual per-capita rate for each
##     source country; Türkiye's hotel occupancy/visitor volume does
##     not affect dengue incidence in Brazil or India.
##
## (2) ×12 scaling factor removed
##     The previous formula multiplied rate_pp_year by 12, which
##     inflated the annual integral by a factor of 12 when S_m is
##     normalized to mean=1 (sum=12). With S_m now in V_in_day,
##     the daily rate is simply rate_pp_year / 365.
##     Verification: Σ_m(daily_rate × days_m) = rate_pp_year/365 × 365
##                 = rate_pp_year ✓
##
## (3) S_m source updated from occupancy to TÜİK departures (v3)
##     See build_seasonality_monthly.R for details.
##
## Thesis formula (Denklem 5):
##   λ_ithal(d,t,m) = V_in(d,t,m) · π_weighted(t,m) · η
##
##   V_in(d,t,m) = V_TR_d × N_total × S_m / 365
##
##   π_weighted(t,m) = 1 - exp(-γ_GBD(t) / 365 × M_climate(t,m) × d_days)
##
##   γ_GBD(t) = Σ_c [ w_c · rate_pp_c(t) ]   (GBD-weighted annual rate)
##
##   References:
##     - Liebig et al. (2019) PLoS ONE — GBD-based formula
##     - Wilder-Smith & Gubler (2008) — travel volume as primary driver
##     - Semenza et al. (2014) PLoS NTD — air travel dengue importation
##     - Stanaway et al. (2016) Lancet Infect Dis — GBD dengue
##     - Massad et al. (2018) Sci Rep — Europe importation model
##
## Inputs (SSP-independent, from DIR_PROCESSED):
##   - population_sentinel_baseline_2024.csv
##   - travel_weights_static.rds
##   - seasonality_monthly.rds  (v3: TÜİK departures)
##
## Inputs (from data_raw/population):
##   - GBD_2023_DATA/GBD_2023_DATA.csv
##   - turizm_verileri/tourist_arrivals_by_country.csv
##
## Inputs (SSP-dependent, from DIR_PROCESSED_SSP):
##   - import_climate_multiplier_monthly.rds
##
## Outputs (SSP-dependent, to DIR_PROCESSED_SSP):
##   - importation_pressure_monthly_2025_2075.csv/.rds
##   - importation_pressure_yearly_2025_2075.csv/.rds
## =========================================================

source("R/01_setup/init.R")

SSP_SCENARIO <- Sys.getenv("SSP_SCENARIO", unset = "ssp245")
cat("\n=============================================\n")
cat("build_importation_pressure_monthly.R (v3 — S_m in V_in)\n")
cat("SSP scenario:", SSP_SCENARIO, "\n")
cat("Started at  :", as.character(Sys.time()), "\n")
cat("=============================================\n\n")

## =========================================================
## 1) Parameters
## =========================================================
BASELINE_YEAR <- 2024L
START_YEAR    <- 2025L
END_YEAR      <- 2075L

ETA    <- 0.25    # Viremic fraction (Wilder-Smith et al. 2019)
D_DAYS <- 3       # Average exposure duration (days)

K_ELASTICITY <- as.numeric(Sys.getenv("K_ELASTICITY",
                                       unset = "0.13"))

GBD_BASELINE_YEARS <- c(2019L, 2021L, 2022L, 2023L)

cat("Parameters:\n")
cat("  ETA =", ETA, "\n")
cat("  D_DAYS =", D_DAYS, "\n")
cat("  K_ELASTICITY =", K_ELASTICITY, "\n")
cat("  GBD baseline years:", paste(GBD_BASELINE_YEARS, collapse = ", "), "\n\n")

norm_key <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- stringr::str_to_lower(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  x
}

## =========================================================
## 2) Load sentinel baseline + travel weights + seasonality
## =========================================================
baseline_path <- file.path(DIR_PROCESSED, "population_sentinel_baseline_2024.csv")
baseline <- readr::read_csv(baseline_path, show_col_types = FALSE)

sentinel <- baseline %>%
  dplyr::filter(year == BASELINE_YEAR) %>%
  dplyr::distinct(district_id, province_name, district_name, .keep_all = TRUE) %>%
  dplyr::select(district_id, province_name, district_name, population) %>%
  dplyr::mutate(
    prov_key = norm_key(province_name),
    dist_key = norm_key(district_name),
    join_key = paste0(prov_key, "_", dist_key)
  )

cat("Sentinel districts:", nrow(sentinel), "\n")

## Travel weights (V_TR per district)
tw_path <- file.path(DIR_PROCESSED, "travel_weights_static.rds")
if (!file.exists(tw_path)) {
  tw_path2 <- file.path(DIR_PROCESSED, "travel_weights.rds")
  if (file.exists(tw_path2)) tw_path <- tw_path2
  else stop("Travel weights file not found.", call. = FALSE)
}
tw <- readRDS(tw_path)

districts <- sentinel %>%
  dplyr::left_join(
    tw %>% dplyr::select(district_id, V_TR = travel_weight),
    by = "district_id"
  )

if (any(is.na(districts$V_TR))) {
  stop("Some districts have no travel weight (V_TR).", call. = FALSE)
}

## Seasonality (v3: TÜİK departures, mean=1)
seas_path <- file.path(DIR_PROCESSED, "seasonality_monthly.rds")
seasonality <- readRDS(seas_path)
cat("Seasonality months:", nrow(seasonality), "\n")
cat("  S_m mean:", round(mean(seasonality$seasonal_multiplier), 4), "\n")
cat("  S_m range:", round(range(seasonality$seasonal_multiplier), 4), "\n")

## =========================================================
## 3) Load M_climate (SSP-dependent)
## =========================================================
mcli_path <- file.path(DIR_PROCESSED_SSP, "import_climate_multiplier_monthly.rds")
USE_M_CLIMATE <- file.exists(mcli_path)

if (USE_M_CLIMATE) {
  mcli <- readRDS(mcli_path)
  cat("M_climate loaded:", nrow(mcli), "rows\n")
  cat("  M_climate range:", round(range(mcli$M_climate), 3), "\n")
} else {
  cat("M_climate not found — using 1.0 for all years.\n")
}

## =========================================================
## 4) Load and process GBD dengue incidence by country
## =========================================================
cat("\n--- Loading GBD 2023 dengue incidence data ---\n")

gbd_path <- file.path(DIR_POP, "GBD_2023_DATA", "GBD_2023_DATA.csv")
if (!file.exists(gbd_path)) {
  stop("GBD data not found at: ", gbd_path,
       "\nExpected: data_raw/population/GBD_2023_DATA/GBD_2023_DATA.csv",
       call. = FALSE)
}

gbd_raw <- readr::read_csv(gbd_path, show_col_types = FALSE)

cat("  GBD raw rows:", nrow(gbd_raw), "\n")
cat("  GBD years:", paste(range(gbd_raw$year), collapse = "–"), "\n")
cat("  GBD locations:", dplyr::n_distinct(gbd_raw$location_name), "\n")

gbd <- gbd_raw %>%
  dplyr::select(
    gbd_country = location_name,
    year,
    incidence_rate = val,
    incidence_upper = upper,
    incidence_lower = lower
  ) %>%
  dplyr::mutate(
    rate_pp       = incidence_rate  / 1e5,
    rate_pp_upper = incidence_upper / 1e5,
    rate_pp_lower = incidence_lower / 1e5
  )

## =========================================================
## 5) Load and process tourist arrivals by country
## =========================================================
cat("\n--- Loading tourist arrivals by country ---\n")

tourist_path <- file.path(DIR_POP, "turizm_verileri",
                          "tourist_arrivals_by_country.csv")
if (!file.exists(tourist_path)) {
  stop("Tourist arrivals data not found at: ", tourist_path, call. = FALSE)
}

tourist_raw <- readr::read_csv(tourist_path, show_col_types = FALSE)

tourist_weights <- tourist_raw %>%
  dplyr::filter(year %in% c(2019, 2021, 2022, 2023),
                !is.na(arrivals), arrivals > 0) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(
    mean_arrivals = mean(arrivals, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    w_c = mean_arrivals / sum(mean_arrivals)
  ) %>%
  dplyr::arrange(dplyr::desc(w_c))

TOTAL_ARRIVALS_ANNUAL <- sum(tourist_weights$mean_arrivals)

cat("  Tourist source countries:", nrow(tourist_weights), "\n")
cat("  Total mean arrivals (annual):", round(TOTAL_ARRIVALS_ANNUAL), "\n")
cat("  Total mean arrivals (daily) :", round(TOTAL_ARRIVALS_ANNUAL / 365), "\n")
cat("  Top 5:\n")
print(head(tourist_weights, 5))

## =========================================================
## 6) Country name harmonization (tourist <-> GBD)
## =========================================================
name_map <- tibble::tribble(
  ~tourist_name,                    ~gbd_name,
  "Russia",                         "Russian Federation",
  "Iran",                           "Iran (Islamic Republic of)",
  "United States",                  "United States of America",
  "Vietnam",                        "Viet Nam",
  "South Korea",                    "Republic of Korea",
  "North Korea",                    "Democratic People's Republic of Korea",
  "Bolivia",                        "Bolivia (Plurinational State of)",
  "Venezuela",                      "Venezuela (Bolivarian Republic of)",
  "Tanzania",                       "United Republic of Tanzania",
  "Democratic Republic of the Congo","Democratic Republic of the Congo",
  "Republic of the Congo",          "Congo",
  "Syria",                          "Syrian Arab Republic",
  "Moldova",                        "Republic of Moldova",
  "Laos",                           "Lao People's Democratic Republic",
  "Ivory Coast",                    "Côte d'Ivoire"
)

tourist_w <- tourist_weights %>%
  dplyr::left_join(name_map, by = c("country" = "tourist_name")) %>%
  dplyr::mutate(
    gbd_match = dplyr::coalesce(gbd_name, country)
  )

matched <- tourist_w %>%
  dplyr::inner_join(
    gbd %>% dplyr::distinct(gbd_country),
    by = c("gbd_match" = "gbd_country")
  )

unmatched <- tourist_w %>%
  dplyr::anti_join(
    gbd %>% dplyr::distinct(gbd_country),
    by = c("gbd_match" = "gbd_country")
  )

cat("\n--- Country matching ---\n")
cat("  Matched:", nrow(matched), "countries (",
    round(sum(matched$w_c) * 100, 1), "% of travel volume)\n")

if (nrow(unmatched) > 0) {
  cat("  Unmatched:", nrow(unmatched), "countries (",
      round(sum(unmatched$w_c) * 100, 1), "% of travel volume)\n")
  cat("  Unmatched top 10:\n")
  print(head(unmatched %>% dplyr::select(country, w_c) %>%
               dplyr::arrange(dplyr::desc(w_c)), 10))
}

matched <- matched %>%
  dplyr::mutate(w_c_norm = w_c / sum(w_c))

cat("  Weights renormalized to matched countries (sum =",
    round(sum(matched$w_c_norm), 4), ")\n")

## =========================================================
## 7) Compute travel-weighted π for GBD years (1990–2023)
## =========================================================
cat("\n--- Computing π_weighted for GBD years ---\n")

gbd_weighted <- matched %>%
  dplyr::select(country, gbd_match, w_c_norm) %>%
  dplyr::inner_join(
    gbd %>% dplyr::select(gbd_country, year, rate_pp, rate_pp_upper, rate_pp_lower),
    by = c("gbd_match" = "gbd_country"),
    relationship = "many-to-many"
  )

pi_annual <- gbd_weighted %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    weighted_rate      = sum(w_c_norm * rate_pp),
    weighted_rate_upper = sum(w_c_norm * rate_pp_upper),
    weighted_rate_lower = sum(w_c_norm * rate_pp_lower),
    n_countries = dplyr::n_distinct(gbd_match),
    .groups = "drop"
  )

cat("  GBD years computed:", nrow(pi_annual), "\n")
cat("  Weighted rate 2019:", signif(pi_annual$weighted_rate[pi_annual$year == 2019], 4), "\n")
cat("  Weighted rate 2023:", signif(pi_annual$weighted_rate[pi_annual$year == 2023], 4), "\n")

## =========================================================
## 8) Project rates beyond GBD (2024–2075) using M_climate
## =========================================================
cat("\n--- Projecting rates for 2024–2075 ---\n")

baseline_rates <- pi_annual %>%
  dplyr::filter(year %in% GBD_BASELINE_YEARS) %>%
  dplyr::summarise(
    base_rate       = mean(weighted_rate),
    base_rate_upper = mean(weighted_rate_upper),
    base_rate_lower = mean(weighted_rate_lower)
  )

cat("  Baseline rate (mean of", paste(GBD_BASELINE_YEARS, collapse = ","), "):",
    signif(baseline_rates$base_rate, 4), "/person/year\n")
cat("  As daily rate:", signif(baseline_rates$base_rate / 365, 4), "/person/day\n")

projection_years <- tibble::tibble(
  year = seq(max(gbd$year) + 1L, END_YEAR, by = 1L)
) %>%
  dplyr::mutate(
    weighted_rate       = baseline_rates$base_rate,
    weighted_rate_upper = baseline_rates$base_rate_upper,
    weighted_rate_lower = baseline_rates$base_rate_lower,
    n_countries = 0L
  )

pi_all_years <- dplyr::bind_rows(
  pi_annual %>% dplyr::filter(year >= START_YEAR),
  projection_years
) %>%
  dplyr::arrange(year)

cat("  Total years:", nrow(pi_all_years),
    "(", sum(pi_all_years$n_countries > 0), "observed,",
    sum(pi_all_years$n_countries == 0), "projected)\n")

pi_scenarios <- dplyr::bind_rows(
  pi_all_years %>%
    dplyr::mutate(pi_scenario = "low",
                  rate_pp_year = weighted_rate_lower),
  pi_all_years %>%
    dplyr::mutate(pi_scenario = "main",
                  rate_pp_year = weighted_rate),
  pi_all_years %>%
    dplyr::mutate(pi_scenario = "high",
                  rate_pp_year = weighted_rate_upper)
) %>%
  dplyr::select(pi_scenario, year, rate_pp_year)

cat("\n  2024 annual rates (per person):\n")
cat("    low :", signif(baseline_rates$base_rate_lower, 4), "\n")
cat("    main:", signif(baseline_rates$base_rate, 4), "\n")
cat("    high:", signif(baseline_rates$base_rate_upper, 4), "\n")

## =========================================================
## 9) Monthly grid
## =========================================================
cat("\n--- Building monthly grid ---\n")

grid <- districts %>%
  dplyr::select(district_id, province_name, district_name, V_TR) %>%
  tidyr::crossing(pi_scenarios %>% dplyr::select(pi_scenario, year, rate_pp_year)) %>%
  tidyr::crossing(seasonality %>% dplyr::select(month, seasonal_multiplier)) %>%
  dplyr::mutate(
    days_in_month = lubridate::days_in_month(
      lubridate::ymd(sprintf("%d-%02d-01", year, month))
    )
  )

# Join M_climate if available
if (USE_M_CLIMATE) {
  grid <- grid %>%
    dplyr::left_join(mcli %>% dplyr::select(year, month, M_climate),
                     by = c("year", "month"))
  if (any(is.na(grid$M_climate))) {
    n_miss <- sum(is.na(grid$M_climate))
    cat("  Filling", n_miss, "missing M_climate values with 1.0\n")
    grid$M_climate[is.na(grid$M_climate)] <- 1.0
  }
} else {
  grid$M_climate <- 1.0
}

## =========================================================
## 10) Compute monthly importation pressure
##
## Formula (Denklem 5, revised):
##   λ_ithal(d,t,m) = V_in(d,m) · π_weighted(t,m) · η
##
## Where:
##   V_in(d,m) = V_TR_d × N_total × S_m / 365
##     → S_m modulates TRAVEL VOLUME (TÜİK departure counts)
##     → V_TR_d is district-level share (TÜİK accommodation data)
##
##   daily_rate = rate_pp_year / 365 × M_climate
##     → rate_pp_year is the GBD-weighted annual per-capita
##       dengue risk across source countries
##     → NO seasonal multiplier here: GBD rate is a country-
##       level annual average; Türkiye's visitor seasonality
##       does not modulate source-country dengue incidence
##
##   π_weighted = 1 - exp(-daily_rate × d_days)
##     → d_days = average exposure duration (3 days)
##
## Annual conservation check:
##   Σ_m(daily_rate × days_m) = rate_pp_year/365 × Σ(days_m)
##                            = rate_pp_year/365 × 365
##                            = rate_pp_year  ✓
##
##   Σ_m(V_in_day_m × days_m) = V_TR × N_total/365 × Σ(S_m × days_m)
##                             ≈ V_TR × N_total/365 × 365   [S_m mean=1]
##                             = V_TR × N_total  ✓
## =========================================================
cat("\n--- Computing importation pressure ---\n")

imp_m <- grid %>%
  dplyr::mutate(
    eta    = ETA,
    d_days = D_DAYS,

    # Daily infection rate for a traveler (GBD-weighted, climate-adjusted)
    # NO S_m here: source-country dengue incidence is independent of
    # Türkiye's seasonal visitor volume
    daily_rate = rate_pp_year / 365 * M_climate,

    # π = probability of infection during d-day exposure (Liebig Eq 1 form)
    pi_weighted = 1 - exp(-daily_rate * d_days),

    # V_in_day: absolute daily arrivals to this district
    # S_m HERE: seasonal variation in international travel volume
    # to Türkiye (TÜİK/EGM departure statistics, mean=1)
    V_in_day = V_TR * TOTAL_ARRIVALS_ANNUAL * seasonal_multiplier / 365,

    # Per-day importation rate
    lambda_import_per_day = V_in_day * pi_weighted * eta,

    # Monthly totals
    lambda_import_window = lambda_import_per_day * pmin(d_days, days_in_month),
    expected_imported_cases_per_month = lambda_import_per_day * days_in_month,

    # Bernoulli: P(at least 1 introduction in the window)
    q_import_month = 1 - exp(-lambda_import_window),

    # For CTMC compatibility
    lambda_import = expected_imported_cases_per_month
  ) %>%
  dplyr::select(
    district_id, province_name, district_name,
    pi_scenario, year, month,
    V_TR, V_in_day, rate_pp_year, daily_rate, pi_weighted, seasonal_multiplier, M_climate,
    eta, days_in_month, d_days,
    lambda_import_per_day,
    lambda_import_window,
    lambda_import,
    expected_imported_cases_per_month,
    q_import_month
  ) %>%
  dplyr::arrange(pi_scenario, district_id, year, month)

## =========================================================
## 11) Diagnostic: expected cases in baseline year
## =========================================================
check_baseline <- imp_m %>%
  dplyr::filter(year == BASELINE_YEAR) %>%
  dplyr::group_by(pi_scenario) %>%
  dplyr::summarise(
    expected_annual_cases = sum(expected_imported_cases_per_month),
    .groups = "drop"
  )

cat("\n--- GBD-weighted expected imported cases (", BASELINE_YEAR, ") ---\n")
print(check_baseline)
cat("  (Note: these are model-derived, not calibration targets)\n")

# Also check 2075 for comparison
check_2075 <- imp_m %>%
  dplyr::filter(year == 2075L) %>%
  dplyr::group_by(pi_scenario) %>%
  dplyr::summarise(
    expected_annual_cases = sum(expected_imported_cases_per_month),
    M_climate_mean = mean(M_climate),
    .groups = "drop"
  )

cat("\n--- Projected imported cases (2075) ---\n")
print(check_2075)

## =========================================================
## 12) Conservation check
## =========================================================
cat("\n--- Annual conservation check (main scenario, 2024) ---\n")
cons_check <- imp_m %>%
  dplyr::filter(pi_scenario == "main", year == BASELINE_YEAR) %>%
  dplyr::group_by(district_id) %>%
  dplyr::summarise(
    sum_daily_rate_x_days = sum(daily_rate * days_in_month),
    rate_pp_year = dplyr::first(rate_pp_year),
    sum_Vin_x_days = sum(V_in_day * days_in_month),
    expected_annual_Vin = dplyr::first(V_TR) * TOTAL_ARRIVALS_ANNUAL,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rate_ratio = sum_daily_rate_x_days / rate_pp_year,
    Vin_ratio = sum_Vin_x_days / expected_annual_Vin
  )

cat("  Rate conservation (should be ~1.0):\n")
print(cons_check %>% dplyr::select(district_id, rate_ratio, Vin_ratio))

## =========================================================
## 13) Save outputs (SSP-specific)
## =========================================================
out_csv <- file.path(DIR_PROCESSED_SSP,
                     "importation_pressure_monthly_2025_2075.csv")
out_rds <- file.path(DIR_PROCESSED_SSP,
                     "importation_pressure_monthly_2025_2075.rds")

readr::write_csv(imp_m, out_csv)
saveRDS(imp_m, out_rds)

# Yearly summary
imp_y <- imp_m %>%
  dplyr::group_by(district_id, province_name, district_name,
                  pi_scenario, year) %>%
  dplyr::summarise(
    lambda_import_year = sum(lambda_import, na.rm = TRUE),
    q_import_year = 1 - prod(1 - q_import_month, na.rm = TRUE),
    M_climate_mean = mean(M_climate, na.rm = TRUE),
    weighted_rate_pp_year = dplyr::first(rate_pp_year),
    .groups = "drop"
  )

out_y_csv <- file.path(DIR_PROCESSED_SSP,
                       "importation_pressure_yearly_2025_2075.csv")
out_y_rds <- file.path(DIR_PROCESSED_SSP,
                       "importation_pressure_yearly_2025_2075.rds")

readr::write_csv(imp_y, out_y_csv)
saveRDS(imp_y, out_y_rds)

# Save diagnostic: country weights and GBD contribution
diag_country <- gbd_weighted %>%
  dplyr::filter(year == 2023L) %>%
  dplyr::mutate(
    contribution_pct = w_c_norm * rate_pp / sum(w_c_norm * rate_pp) * 100
  ) %>%
  dplyr::select(country, gbd_match, w_c_norm, rate_pp,
                contribution_pct) %>%
  dplyr::arrange(dplyr::desc(contribution_pct))

diag_path <- file.path(DIR_PROCESSED_SSP,
                       "importation_country_contributions.csv")
readr::write_csv(diag_country, diag_path)
cat("\n  Country contributions saved:", diag_path, "\n")

cat("\n=============================================\n")
cat("Finished at   :", as.character(Sys.time()), "\n")
cat("SSP scenario  :", SSP_SCENARIO, "\n")
cat("M_climate used:", USE_M_CLIMATE, "\n")
cat("GBD data      : country-specific, under-report corrected\n")
cat("ρ             : NOT USED (GBD already corrected)\n")
cat("S_m location  : V_in(t) (travel volume, NOT γ_GBD)\n")
cat("S_m source    : TÜİK departures (v3)\n")
cat("Output CSV    :", out_csv, "\n")
cat("Output RDS    :", out_rds, "\n")
cat("Rows          :", nrow(imp_m), "\n")
cat("Districts     :", dplyr::n_distinct(imp_m$district_id), "\n")
cat("Years         :", min(imp_m$year), "-", max(imp_m$year), "\n")
cat("Scenarios     :", paste(unique(imp_m$pi_scenario), collapse = ", "), "\n")
cat("=============================================\n")

message("\n=== DONE (", SSP_SCENARIO, ") ===")
