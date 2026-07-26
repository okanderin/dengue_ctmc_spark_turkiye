## =========================================================
## Project : r_project_tez
## Script  : R/02_data/build_importation_pressure_monthly.R
## Purpose : Monthly importation pressure dataset (sentinel districts)
##           SSP-aware: M_climate varies by scenario.
##
## VERSION 4 — METHOD 1: district-level RAW arrivals (A_d), direct
##
## KEY CHANGE FROM v3
##   V_in_day now uses each district's OWN raw foreign accommodation
##   arrivals A_d directly, instead of allocating the NATIONAL total
##   across the 5 sentinels via a within-5-district share (V_TR).
##
##     v3 (K1):  V_in_day = V_TR_d × N_total × S_m / 365
##               V_TR_d = A_d / Σ_{5} A_j   ⟹  Σ_{5} V_in = N_total
##               (all Türkiye funneled into 5 districts — dimensional error)
##
##     v4 (M1):  V_in_day = A_d × S_m / 365
##               A_d = TÜİK district foreign accommodation arrivals
##               (each district gets its OWN arrivals; no funneling)
##
##   All other terms (π_weighted, η, λ_import, downstream P_ufuk) are
##   UNCHANGED. v4 = v3 × c, with c = (Σ_{5} A_j)/N_total a single global
##   scalar ⟹ district ranking, threshold years and season length are
##   INVARIANT; only the absolute P_ufuk scale changes.
##   (Proof + numerical check: invariance_check_importation.R)
##
##   Source-country composition (w_c) remains NATIONAL — an explicit
##   assumption that the source mix is similar across sentinels.
##
## Thesis formula (Denklem 5, revised):
##   λ_ithal(d,t,m) = V_in(d,m) · π_weighted(t,m) · η
##   V_in(d,m)      = A_d × S_m / 365
##   π_weighted     = 1 - exp(-γ_GBD(t,m) · d_days)
##   γ_GBD(t,m)     = (GBD_insidans(t)/100000) · M_climate(t,m)
##
## References:
##   - Liebig et al. (2019) PLoS ONE  [115]  GBD-based arrival model
##   - Semenza et al. (2014) PLoS NTD [104]  air-travel importation
##   - Massad et al. (2018) Sci Rep   [113]  introduction probability
##   - Wilder-Smith & Gubler (2008)   [116]  travel as primary driver
##
## Inputs (SSP-independent, from DIR_PROCESSED):
##   - population_sentinel_baseline_2024.csv
##   - travel_weights_static.rds        (MUST contain `gelis_yabanci` = A_d)
##   - seasonality_monthly.rds          (v3: TÜİK departures)
## Inputs (from data_raw/population):
##   - GBD_2023_DATA/GBD_2023_DATA.csv
##   - turizm_verileri/tourist_arrivals_by_country.csv
## Inputs (SSP-dependent, from DIR_PROCESSED_SSP):
##   - import_climate_multiplier_monthly.rds
## Outputs (SSP-dependent, to DIR_PROCESSED_SSP):
##   - importation_pressure_monthly_2025_2075.csv/.rds
##   - importation_pressure_yearly_2025_2075.csv/.rds
## =========================================================

source("R/01_setup/init.R")

SSP_SCENARIO <- Sys.getenv("SSP_SCENARIO", unset = "ssp245")
cat("\n=============================================\n")
cat("build_importation_pressure_monthly.R (v4 — METHOD 1: A_d direct)\n")
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

## Travel weights + RAW arrivals A_d per district
## -------------------------------------------------------------------
## METHOD 1: pull the raw `gelis_yabanci` column (A_d) directly.
## `travel_weight` (V_TR) is retained for reference/diagnostics only.
tw_path <- file.path(DIR_PROCESSED, "travel_weights_zone.rds")
if (!file.exists(tw_path)) {
  tw_path2 <- file.path(DIR_PROCESSED, "travel_weights_zone.rds")
  if (file.exists(tw_path2)) tw_path <- tw_path2
  else stop("Travel weights file not found.", call. = FALSE)
}
tw <- readRDS(tw_path)

if (!"gelis_yabanci" %in% names(tw)) {
  stop("travel_weights_zone.rds `gelis_yabanci` (A_d) sutununu icermiyor; ",
       "build_travel_weights.R'yi (v2) yeniden calistirin.", call. = FALSE)
}

districts <- sentinel %>%
  dplyr::left_join(
    tw %>% dplyr::select(
      district_id,
      V_TR = travel_weight,     # korunuyor: goreli 5-ilce payi (yalnizca tanilama)
      A_d  = gelis_yabanci      # METHOD 1: ham yillik yabanci tesis gelisi
    ),
    by = "district_id"
  )

if (any(is.na(districts$V_TR)) || any(is.na(districts$A_d))) {
  stop("Bazi ilcelerde V_TR ya da A_d (gelis_yabanci) eksik.", call. = FALSE)
}

cat("  A_d (ham yabanci gelis) yuklendi:\n")
print(districts %>% dplyr::select(district_name, A_d, V_TR))

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
##    (national source-country weights w_c; N_total for diagnostics)
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

## NOTE (METHOD 1): TOTAL_ARRIVALS_ANNUAL is NO LONGER used in V_in.
## Kept for source-country diagnostics and for computing the invariance
## scale c = sum(A_d)/N_total (see invariance_check_importation.R).
TOTAL_ARRIVALS_ANNUAL <- sum(tourist_weights$mean_arrivals)

cat("  Tourist source countries:", nrow(tourist_weights), "\n")
cat("  Total mean arrivals (annual):", round(TOTAL_ARRIVALS_ANNUAL), "\n")
cat("  Total mean arrivals (daily) :", round(TOTAL_ARRIVALS_ANNUAL / 365), "\n")
cat("  Top 5:\n")
print(head(tourist_weights, 5))

## Invariance scale (informational): NEW = OLD × c, with c constant.
c_scale <- sum(districts$A_d) / TOTAL_ARRIVALS_ANNUAL
cat(sprintf("\n  [METHOD 1] c = sum(A_d)/N_total = %.5f  (v3 Lambda ~%.1f kat sisikti)\n",
            c_scale, 1 / c_scale))

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
  dplyr::select(district_id, province_name, district_name, V_TR, A_d) %>%
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
## Formula (Denklem 5, METHOD 1):
##   λ_ithal(d,t,m) = V_in(d,m) · π_weighted(t,m) · η
##
## Where:
##   V_in(d,m) = A_d × S_m / 365      (Yontem 1: dogrudan ilce gelisi)
##     → A_d : TUIK ilce bazli yillik yabanci tesis gelisi (ham)
##     → S_m : uluslararasi seyahat hacminin aylik degisimi (ortalama=1)
##     → NO N_total, NO 5-ilce normalizasyonu
##
##   daily_rate = rate_pp_year / 365 × M_climate
##     → rate_pp_year is the GBD-weighted annual per-capita dengue risk
##       across NATIONAL source countries
##     → NO S_m here: source-country dengue incidence is independent of
##       Türkiye's visitor seasonality
##
##   π_weighted = 1 - exp(-daily_rate × d_days)    (d_days = 3)
##
## Annual conservation check:
##   Σ_m(daily_rate × days_m) = rate_pp_year/365 × 365 = rate_pp_year  ✓
##   Σ_m(V_in_day_m × days_m) = A_d/365 × Σ(S_m × days_m)
##                            ≈ A_d/365 × 365 = A_d                     ✓
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

    # V_in_day: bu ilcenin MUTLAK gunluk yabanci gelisi (METHOD 1)
    # Dogrudan ham ilce verisi A_d; ulusal toplam / normalizasyon YOK.
    # S_m: uluslararasi seyahat hacminin aylik mevsimselligi (ortalama=1)
    V_in_day = A_d * seasonal_multiplier / 365,

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
    V_TR, A_d, V_in_day, rate_pp_year, daily_rate, pi_weighted, seasonal_multiplier, M_climate,
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
## 12) Conservation check (METHOD 1: annual V_in must equal A_d)
## =========================================================
cat("\n--- Annual conservation check (main scenario, 2024) ---\n")
cons_check <- imp_m %>%
  dplyr::filter(pi_scenario == "main", year == BASELINE_YEAR) %>%
  dplyr::group_by(district_id) %>%
  dplyr::summarise(
    sum_daily_rate_x_days = sum(daily_rate * days_in_month),
    rate_pp_year          = dplyr::first(rate_pp_year),
    sum_Vin_x_days        = sum(V_in_day * days_in_month),
    expected_annual_Vin   = dplyr::first(A_d),   # METHOD 1: yillik V_in = A_d
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rate_ratio = sum_daily_rate_x_days / rate_pp_year,
    Vin_ratio  = sum_Vin_x_days / expected_annual_Vin   # ~1.0 olmali
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
cat("V_in method   : METHOD 1 — raw district arrivals A_d (no N_total)\n")
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
