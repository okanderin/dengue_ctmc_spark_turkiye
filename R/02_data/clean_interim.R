## =========================================================
## Project : r_project_tez
## Script  : R/02_data/clean_interim.R
## Purpose : Clean interim climate (district-month) dataset for modeling
## Input   : DIR_INTERIM_SSP/climate_district_monthly_2015_2100.rds
## Output  : DIR_PROCESSED_SSP/climate_sentinel_monthly_long_2015_2100.rds
##         : DIR_PROCESSED_SSP/climate_sentinel_monthly_wide_2015_2100.rds
## Log     : DIR_LOGS/clean_interim_{ssp}_log.txt
## =========================================================

## ---- 0) Initialize (packages + paths + options)
source("R/01_setup/init.R")

## ---------------------------------------------------------
## 1) Logging
## ---------------------------------------------------------
log_file <- file.path(DIR_LOGS, paste0("clean_interim_", SSP_SCENARIO, "_log.txt"))
sink(log_file, split = TRUE)
on.exit({ try(sink(), silent = TRUE) }, add = TRUE)

cat("=============================================\n")
cat("clean_interim.R started at :", as.character(Sys.time()), "\n")
cat("Project ROOT               :", ROOT, "\n")
cat("SSP_SCENARIO               :", SSP_SCENARIO, "\n")
cat("DIR_INTERIM_SSP            :", DIR_INTERIM_SSP, "\n")
cat("DIR_PROCESSED_SSP          :", DIR_PROCESSED_SSP, "\n")
cat("=============================================\n\n")

## ---------------------------------------------------------
## 2) Read interim dataset
## ---------------------------------------------------------
in_file <- fp_climate_interim
if (!file.exists(in_file)) stop("Missing input: ", in_file, call. = FALSE)

x <- readRDS(in_file)

cat("Loaded:", in_file, "\n")
cat("Rows:", nrow(x), "Cols:", ncol(x), "\n")
cat("Columns:", paste(names(x), collapse = ", "), "\n\n")

## ---------------------------------------------------------
## 3) Basic validation
## ---------------------------------------------------------
required_cols <- c("district_id", "province_name", "district_name",
                   "year", "month", "date", "variable", "value")
miss_cols <- setdiff(required_cols, names(x))
if (length(miss_cols) > 0) {
  stop("Missing columns: ", paste(miss_cols, collapse = ", "), call. = FALSE)
}

# Expected variables — uses VAR_HUMIDITY from paths.R (= "hurs")
vars_expected <- c("tas", "tasmin", "tasmax", "pr", VAR_HUMIDITY)
vars_have <- sort(unique(x$variable))
cat("Expected variables:", paste(vars_expected, collapse = ", "), "\n")
cat("Found variables   :", paste(vars_have, collapse = ", "), "\n")

if (!all(vars_expected %in% vars_have)) {
  missing_vars <- setdiff(vars_expected, vars_have)
  stop("Missing expected variables: ", paste(missing_vars, collapse = ", "),
       "\nAvailable: ", paste(vars_have, collapse = ", "), call. = FALSE)
}

## ---------------------------------------------------------
## 4) Sentinel table (hard-coded, already validated)
## ---------------------------------------------------------
sentinel <- tibble::tibble(
  prov = c("Artvin", "Zinguldak", "Istanbul", "Isparta", "Mugla"),
  dist = c("Hopa",   "Merkez",    "Kartal",  "Eğirdir", "Fethiye")
) %>%
  dplyr::mutate(
    prov_norm = stringr::str_to_lower(
      stringi::stri_trans_general(prov, "Latin-ASCII")),
    dist_norm = stringr::str_to_lower(
      stringi::stri_trans_general(dist, "Latin-ASCII"))
  )

## Normalize names in dataset (to join robustly)
x2 <- x %>%
  dplyr::mutate(
    prov_norm = stringr::str_to_lower(
      stringi::stri_trans_general(province_name, "Latin-ASCII")),
    dist_norm = stringr::str_to_lower(
      stringi::stri_trans_general(district_name, "Latin-ASCII"))
  )

## Sentinel join check
have <- x2 %>% dplyr::distinct(prov_norm, dist_norm)
anti <- sentinel %>%
  dplyr::anti_join(have, by = c("prov_norm", "dist_norm"))

if (nrow(anti) > 0) {
  stop("Sentinel mismatch (not found in dataset):\n",
       paste(paste0(anti$prov, " / ", anti$dist), collapse = "\n"),
       call. = FALSE)
}

## Filter to sentinel districts (should stay 5)
x_s <- x2 %>%
  dplyr::semi_join(sentinel, by = c("prov_norm", "dist_norm")) %>%
  dplyr::select(-prov_norm, -dist_norm)

cat("Sentinel districts retained:", dplyr::n_distinct(x_s$district_id), "\n")

## ---------------------------------------------------------
## 4b) Standardize variable names for downstream compatibility
##     CNRM-CM6-1-HR uses "hurs" (near-surface relative humidity).
##     Model scripts (ctmc_spark.R, parameter_functions.R) expect "hur".
##     Rename here so all downstream code works without changes.
## ---------------------------------------------------------
if ("hurs" %in% unique(x_s$variable)) {
  x_s <- x_s %>% dplyr::mutate(variable = ifelse(variable == "hurs", "hur", variable))
  cat("Renamed variable: hurs -> hur (downstream compatibility)\n")
}

## Update vars_expected to match the standardized names
vars_expected <- c("tas", "tasmin", "tasmax", "pr", "hur")

## ---------------------------------------------------------
## 5) Standardize monthly date key
##    NetCDF monthly times often come as mid-month (e.g., 16th).
##    For modeling, set to first day of month.
## ---------------------------------------------------------
x_s <- x_s %>%
  dplyr::mutate(
    ym   = sprintf("%04d-%02d", year, month),
    date = as.Date(paste0(ym, "-01"))
  )

cat("Date range (standardized):",
    as.character(min(x_s$date)), "->",
    as.character(max(x_s$date)), "\n")

## ---------------------------------------------------------
## 6) Consistency checks (months per district-variable)
## ---------------------------------------------------------
expected_months <- length(
  seq(as.Date("2015-01-01"), as.Date("2100-12-01"), by = "1 month")
)

miss_check <- x_s %>%
  dplyr::count(district_id, variable, name = "n_obs") %>%
  dplyr::mutate(missing = expected_months - n_obs) %>%
  dplyr::filter(missing != 0)

if (nrow(miss_check) > 0) {
  cat("\nWARNING: Missing months detected:\n")
  print(miss_check, n = Inf)
  stop("Missing months exist; investigate NetCDF coverage.", call. = FALSE)
} else {
  cat("OK: No missing months per district-variable.\n")
}

## ---------------------------------------------------------
## 7) Apply bias correction (Delta method)
##    Correction factors are pre-computed by build_bias_correction.R
##    and stored in data_processed/bias_correction_factors.rds
## ---------------------------------------------------------
fp_bc <- file.path(DIR_PROCESSED_SHARED, "bias_correction_factors.rds")

if (file.exists(fp_bc)) {
  cat("\n>>> Applying bias correction (Delta method) ...\n")

  bc <- readRDS(fp_bc) %>%
    dplyr::select(district_id, month, variable, method, correction)

  # Map correction factor variable names to current data variable names
  # In correction factors: "hurs"; in data after step 4b: "hur"
  bc <- bc %>%
    dplyr::mutate(variable = dplyr::if_else(variable == "hurs", "hur", variable))

  n_before <- nrow(x_s)

  x_s <- x_s %>%
    dplyr::left_join(bc, by = c("district_id", "month", "variable"))

  # Apply correction
  x_s <- x_s %>%
    dplyr::mutate(
      value_raw = value,
      value = dplyr::case_when(
        is.na(correction)        ~ value,                 # no correction available
        method == "additive"     ~ value + correction,
        method == "multiplicative" ~ value * correction,
        TRUE                     ~ value
      ),
      # Bound humidity to [0, 100] after additive correction
      value = dplyr::if_else(variable == "hur",
                              pmin(pmax(value, 0), 100),
                              value),
      # Bound precipitation to >= 0 after multiplicative correction
      value = dplyr::if_else(variable == "pr",
                              pmax(value, 0),
                              value)
    )

  # Report correction magnitude
  bc_summary <- x_s %>%
    dplyr::filter(!is.na(correction)) %>%
    dplyr::group_by(variable, method) %>%
    dplyr::summarise(
      mean_raw       = round(mean(value_raw, na.rm = TRUE), 2),
      mean_corrected = round(mean(value, na.rm = TRUE), 2),
      mean_shift     = round(mean(value - value_raw, na.rm = TRUE), 3),
      .groups = "drop"
    )

  cat("  Bias correction applied.\n")
  print(bc_summary)

  x_s <- x_s %>% dplyr::select(-value_raw, -method, -correction)

  stopifnot(nrow(x_s) == n_before)
  cat("  Row count preserved:", nrow(x_s), "\n")
} else {
  cat("\n>>> SKIP bias correction: ", fp_bc, " not found.\n")
  cat("  Run build_bias_correction.R first if correction is desired.\n")
}

## ---------------------------------------------------------
## 8) Create outputs
##    a) Long: district-month-variable
##    b) Wide: one row per district-month, variables as columns
## ---------------------------------------------------------
out_long <- x_s %>%
  dplyr::arrange(district_id, date, variable)

out_wide <- out_long %>%
  dplyr::select(district_id, province_name, district_name,
                year, month, date, variable, value) %>%
  tidyr::pivot_wider(
    names_from = variable,
    values_from = value
  ) %>%
  dplyr::arrange(district_id, date)

## Quick sanity on wide
stopifnot(all(vars_expected %in% names(out_wide)))

## ---------------------------------------------------------
## 9) Save processed outputs (SSP-specific directory)
## ---------------------------------------------------------
if (!dir.exists(DIR_PROCESSED_SSP)) {
  dir.create(DIR_PROCESSED_SSP, recursive = TRUE)
}

out_long_file <- file.path(DIR_PROCESSED_SSP,
                           "climate_sentinel_monthly_long_2015_2100.rds")
out_wide_file <- file.path(DIR_PROCESSED_SSP,
                           "climate_sentinel_monthly_wide_2015_2100.rds")

saveRDS(out_long, out_long_file)
saveRDS(out_wide, out_wide_file)

cat("\n=============================================\n")
cat("clean_interim.R finished at :", as.character(Sys.time()), "\n")
cat("SSP scenario               :", SSP_SCENARIO, "\n")
cat("Saved long :", out_long_file, "\n")
cat("Saved wide :", out_wide_file, "\n")
cat("Long rows  :", nrow(out_long), "\n")
cat("Wide rows  :", nrow(out_wide), "\n")
cat("Districts  :", dplyr::n_distinct(out_wide$district_id), "\n")
cat("=============================================\n")

message("\n=== DONE (", SSP_SCENARIO, ") ===")
message("Saved long: ", out_long_file)
message("Saved wide: ", out_wide_file)
