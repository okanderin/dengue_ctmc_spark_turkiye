## =========================================================
## Project : r_project_tez
## Script  : R/02_data/build_bias_correction.R
## Purpose : Compute bias correction factors for CMIP6
##           projections using ERA5-Land reanalysis as reference.
##
## Method  : Monthly Delta Method (per district, per month-of-year)
##           - Temperature (tas): additive
##           - Precipitation (pr): multiplicative
##           - Humidity (hurs): additive, bounded [0,100]
##
## Note    : ERA5-Land only provides t2m, d2m, tp — no tasmin/tasmax.
##           Therefore bias correction is applied ONLY to:
##             tas  (CMIP6) vs t2m (ERA5)
##             pr   (CMIP6) vs tp  (ERA5)
##             hurs (CMIP6) vs RH computed from t2m & d2m (ERA5)
##           tasmin and tasmax pass through uncorrected.
##
## Reference: Bagcaci et al. (2021) doi:10.1016/j.atmosres.2021.105576
##            — validated Delta method for CMIP6 over Turkey with ERA ref
##
## Inputs  :
##   DIR_ERA5          : era5_land_1981_2014.nc  (CDS monthly means)
##   DIR_CLIMATE_HIST  : cmip6_hist_1981_2014/{hurs,pr,tas}/*.nc
##   DIR_SHP           : gadm41_TUR_2.shp
##
## Outputs :
##   data_processed/bias_correction_factors.rds
##   data_processed/bias_correction_factors.csv
##   data_processed/bias_correction_diagnostics.csv
##
## Usage   : Source ONCE before running import_raw.R
## =========================================================

source("R/01_setup/init.R")

## ---------------------------------------------------------
## 0) Logging
## ---------------------------------------------------------
log_file <- file.path(DIR_LOGS, "build_bias_correction_log.txt")
sink(log_file, split = TRUE)
on.exit({ try(sink(), silent = TRUE) }, add = TRUE)

cat("=============================================\n")
cat("build_bias_correction.R started at:", as.character(Sys.time()), "\n")
cat("DIR_ERA5         :", DIR_ERA5, "\n")
cat("DIR_CLIMATE_HIST :", DIR_CLIMATE_HIST, "\n")
cat("DIR_SHP          :", DIR_SHP, "\n")
cat("=============================================\n\n")

## ---------------------------------------------------------
## 1) Load sentinel district polygons
## ---------------------------------------------------------
shp_path <- file.path(DIR_SHP, "gadm41_TUR_2.shp")
stopifnot(file.exists(shp_path))

districts_raw <- sf::st_read(shp_path, quiet = TRUE)

geom_col <- attr(districts_raw, "sf_column")
if (!is.null(geom_col) && geom_col != "geometry") {
  names(districts_raw)[names(districts_raw) == geom_col] <- "geometry"
  attr(districts_raw, "sf_column") <- "geometry"
}

norm_tr <- function(x) {
  stringi::stri_trans_general(stringr::str_to_lower(stringr::str_squish(x)),
                              "Latin-ASCII")
}

sentinel <- tibble::tibble(
  prov = c("Artvin", "Zinguldak", "Istanbul", "Isparta", "Mugla"),
  dist = c("Hopa",   "Merkez",    "Kartal",   "Egirdir", "Fethiye")
) %>%
  dplyr::mutate(prov_norm = norm_tr(prov), dist_norm = norm_tr(dist))

districts <- districts_raw %>%
  dplyr::mutate(prov_norm = norm_tr(NAME_1), dist_norm = norm_tr(NAME_2)) %>%
  dplyr::inner_join(sentinel %>% dplyr::select(prov_norm, dist_norm),
                    by = c("prov_norm", "dist_norm")) %>%
  dplyr::mutate(
    province_name = NAME_1,
    district_name = NAME_2,
    district_id   = as.character(GID_2)
  ) %>%
  dplyr::select(district_id, province_name, district_name, geometry) %>%
  sf::st_make_valid() %>%
  sf::st_transform(4326)

stopifnot(nrow(districts) == 5)
cat("Sentinel districts loaded:", nrow(districts), "\n\n")

## ---------------------------------------------------------
## 2) Helper: parse dates from raster (handles ERA5-Land)
##    ERA5-Land CDS layer names: "varname_valid_time=EPOCH_SECONDS"
## ---------------------------------------------------------
get_raster_dates <- function(r, fallback_start = "1981-01-01") {
  n <- terra::nlyr(r)

  # Try 1: terra::time()
  tt <- tryCatch(terra::time(r), error = function(e) NULL)
  if (!is.null(tt) && !all(is.na(tt))) {
    return(as.Date(tt))
  }

  # Try 2: Unix epoch in layer names
  epochs <- as.numeric(stringr::str_extract(names(r), "\\d+$"))
  if (!all(is.na(epochs))) {
    dates <- as.Date(as.POSIXct(epochs, origin = "1970-01-01", tz = "UTC"))
    if (length(dates) == n) {
      cat("    [time] parsed from epoch in layer names\n")
      return(dates)
    }
  }

  # Try 3: synthetic monthly sequence
  cat("    [time] WARNING: generating synthetic monthly dates from",
      fallback_start, "\n")
  seq(as.Date(fallback_start), by = "month", length.out = n)
}

## ---------------------------------------------------------
## 3) Helper: extract district means from a raster
## ---------------------------------------------------------
extract_district_monthly <- function(r, districts_sf, varname) {

  n_layers <- terra::nlyr(r)
  dates    <- get_raster_dates(r)
  stopifnot(length(dates) == n_layers)

  # CRS
  if (!grepl("EPSG:4326|WGS 84|longlat", terra::crs(r), ignore.case = TRUE)) {
    r <- terra::project(r, "EPSG:4326")
  }

  # Unit conversion
  if (varname %in% c("tas", "tasmin", "tasmax", "t2m", "d2m")) {
    med <- terra::global(r[[1]], "mean", na.rm = TRUE)[[1]]
    if (!is.na(med) && med > 100) {
      r <- r - 273.15
      cat("    [units] ", varname, " Kelvin -> Celsius\n")
    }
  } else if (varname %in% c("pr", "tp")) {
    med <- terra::global(r[[1]], "mean", na.rm = TRUE)[[1]]
    if (!is.na(med) && med < 0.1) {
      # ERA5-Land monthly tp: metres of water -> mm
      r <- r * 1000
      cat("    [units] ", varname, " metres -> mm (monthly total)\n")
    } else if (!is.na(med) && med < 1) {
      # CMIP6 pr: kg/m2/s -> mm/day
      r <- r * 86400
      cat("    [units] ", varname, " kg/m2/s -> mm/day\n")
    }
  }

  # Extract
  m <- exactextractr::exact_extract(r, districts_sf, fun = "mean", progress = FALSE)

  # Build long table
  # exact_extract returns:
  #   multi-layer: data.frame (n_districts x n_layers)
  #   single-layer: numeric vector (n_districts)
  if (is.numeric(m) && !is.data.frame(m)) {
    m_df <- data.frame(mean = m)
  } else {
    m_df <- as.data.frame(m)
  }
  n_dist <- nrow(m_df)
  n_cols <- ncol(m_df)

  # Safety: number of columns must match number of layers
  if (n_cols != n_layers) {
    warning("exact_extract returned ", n_cols, " columns but expected ",
            n_layers, " layers for ", varname, ". Attempting fallback.")
    # Fallback: terra::extract
    m2 <- terra::extract(r, terra::vect(districts_sf), fun = mean, na.rm = TRUE)
    m_df <- as.data.frame(m2[, -1])  # drop ID column
    n_cols <- ncol(m_df)
  }

  result <- list()
  for (d in seq_len(n_dist)) {
    vals <- as.numeric(m_df[d, ])
    result[[d]] <- tibble::tibble(
      district_id = districts_sf$district_id[d],
      date  = dates[seq_len(n_cols)],
      year  = lubridate::year(dates[seq_len(n_cols)]),
      month = lubridate::month(dates[seq_len(n_cols)]),
      value = vals
    )
  }

  dplyr::bind_rows(result)
}

## ---------------------------------------------------------
## 4) Extract ERA5-Land data (1981-2014)
##    ERA5-Land CDS multi-variable NetCDF: varnames() returns
##    one entry per variable (not per layer). We split layers
##    by their name prefix: "d2m_valid_time=...", "t2m_...", "tp_..."
## ---------------------------------------------------------
cat(">>> Extracting ERA5-Land (1981-2014) ...\n")

era5_nc <- list.files(DIR_ERA5, pattern = "\\.nc$", full.names = TRUE)
stopifnot(length(era5_nc) > 0)

r_era5 <- terra::rast(era5_nc[1])
lyr_names <- names(r_era5)

cat("  File:", basename(era5_nc[1]), "\n")
cat("  Total layers:", terra::nlyr(r_era5), "\n")

# Extract variable prefix from layer names (e.g. "d2m" from "d2m_valid_time=347155200")
lyr_prefix <- stringr::str_extract(lyr_names, "^[^_]+")
unique_prefixes <- unique(lyr_prefix)
cat("  Detected layer prefixes:", paste(unique_prefixes, collapse = ", "), "\n")
cat("  Layers per prefix:\n")
print(table(lyr_prefix))
cat("\n")

# Split by prefix
era5_by_var <- list()
for (vn in unique_prefixes) {
  idx <- which(lyr_prefix == vn)
  r_sub <- r_era5[[idx]]
  cat("  Extracting '", vn, "' (", length(idx), " layers) ...\n", sep = "")
  era5_by_var[[vn]] <- extract_district_monthly(r_sub, districts, vn)
}

# Compute RH from t2m and d2m (Magnus formula)
cat("\n  Computing RH from t2m and d2m ...\n")
stopifnot("t2m" %in% names(era5_by_var), "d2m" %in% names(era5_by_var))

ta_df <- era5_by_var[["t2m"]] %>%
  dplyr::select(district_id, year, month, ta = value)
dp_df <- era5_by_var[["d2m"]] %>%
  dplyr::select(district_id, year, month, dp = value)

rh_df <- dplyr::inner_join(ta_df, dp_df,
                           by = c("district_id", "year", "month")) %>%
  dplyr::mutate(
    a  = 17.625,
    b  = 243.04,
    rh = 100 * exp((a * dp) / (b + dp)) / exp((a * ta) / (b + ta)),
    rh = pmin(pmax(rh, 0), 100)
  )

cat("  RH range:", round(min(rh_df$rh, na.rm = TRUE), 1), " - ",
    round(max(rh_df$rh, na.rm = TRUE), 1), "%\n\n")

# Canonical ERA5 table
era5 <- dplyr::bind_rows(
  era5_by_var[["t2m"]] %>% dplyr::mutate(variable = "tas"),
  rh_df %>%
    dplyr::transmute(district_id, year, month, value = rh, variable = "hurs"),
  era5_by_var[["tp"]] %>% dplyr::mutate(variable = "pr")
) %>%
  dplyr::filter(year >= 1981, year <= 2014)

cat("  ERA5 canonical table:", nrow(era5), "rows |",
    paste(unique(era5$variable), collapse = ", "), "\n\n")

## ---------------------------------------------------------
## 5) Extract CMIP6-Historical data (1981-2014)
## ---------------------------------------------------------
cat(">>> Extracting CMIP6-Historical (1981-2014) ...\n")

hist_vars <- c(tas = "tas", pr = "pr", hurs = "hurs")

hist_all <- list()
for (i in seq_along(hist_vars)) {
  vname  <- names(hist_vars)[i]
  folder <- hist_vars[i]
  v_dir  <- file.path(DIR_CLIMATE_HIST, folder)

  if (!dir.exists(v_dir)) {
    cat("  SKIP (not found):", v_dir, "\n")
    next
  }

  nc_files <- list.files(v_dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) {
    cat("  SKIP (no .nc):", v_dir, "\n")
    next
  }

  cat("  Extracting '", vname, "' from ", length(nc_files), " file(s) ...\n", sep = "")
  r <- terra::rast(nc_files)
  hist_all[[vname]] <- extract_district_monthly(r, districts, vname) %>%
    dplyr::mutate(variable = vname)
}

cmip6_hist <- dplyr::bind_rows(hist_all) %>%
  dplyr::filter(year >= 1981, year <= 2014)

cat("  CMIP6-hist:", nrow(cmip6_hist), "rows |",
    paste(unique(cmip6_hist$variable), collapse = ", "), "\n\n")

## ---------------------------------------------------------
## 6) Compute monthly climatologies
## ---------------------------------------------------------
cat(">>> Computing climatologies (1981-2014) ...\n")

clim_era5 <- era5 %>%
  dplyr::group_by(district_id, month, variable) %>%
  dplyr::summarise(era5_mean = mean(value, na.rm = TRUE),
                   era5_sd   = sd(value, na.rm = TRUE),
                   n_era5    = dplyr::n(), .groups = "drop")

clim_hist <- cmip6_hist %>%
  dplyr::group_by(district_id, month, variable) %>%
  dplyr::summarise(hist_mean = mean(value, na.rm = TRUE),
                   hist_sd   = sd(value, na.rm = TRUE),
                   n_hist    = dplyr::n(), .groups = "drop")

## ---------------------------------------------------------
## 7) Compute correction factors
## ---------------------------------------------------------
cat(">>> Computing correction factors ...\n")

correction_type <- c(tas = "additive", pr = "multiplicative", hurs = "additive")

bc_factors <- dplyr::inner_join(clim_era5, clim_hist,
                                by = c("district_id", "month", "variable")) %>%
  dplyr::mutate(
    method = correction_type[variable],
    delta_add  = era5_mean - hist_mean,
    ratio_mult = dplyr::if_else(hist_mean > 0.001, era5_mean / hist_mean, 1.0),
    correction = dplyr::case_when(
      method == "additive"       ~ delta_add,
      method == "multiplicative" ~ ratio_mult,
      TRUE                       ~ 0
    )
  )

cat("  Factors computed:", nrow(bc_factors), "rows\n")

# Summary
cat("\n--- Bias Correction Summary ---\n")
bc_factors %>%
  dplyr::group_by(variable, method) %>%
  dplyr::summarise(
    mean_corr      = round(mean(correction, na.rm = TRUE), 3),
    min_corr       = round(min(correction, na.rm = TRUE), 3),
    max_corr       = round(max(correction, na.rm = TRUE), 3),
    mean_era5      = round(mean(era5_mean, na.rm = TRUE), 2),
    mean_cmip6hist = round(mean(hist_mean, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  print(n = Inf)

cat("\n  Note: tasmin, tasmax pass through UNCORRECTED\n")

## ---------------------------------------------------------
## 8) Save
## ---------------------------------------------------------
out_rds  <- file.path(DIR_PROCESSED_SHARED, "bias_correction_factors.rds")
out_csv  <- file.path(DIR_PROCESSED_SHARED, "bias_correction_factors.csv")
diag_csv <- file.path(DIR_PROCESSED_SHARED, "bias_correction_diagnostics.csv")

saveRDS(bc_factors, out_rds)
readr::write_csv(bc_factors, out_csv)

diag_df <- dplyr::full_join(
  era5 %>% dplyr::rename(era5_value = value),
  cmip6_hist %>% dplyr::rename(hist_value = value),
  by = c("district_id", "year", "month", "variable")
) %>%
  dplyr::mutate(bias = era5_value - hist_value)
readr::write_csv(diag_df, diag_csv)

cat("\n=============================================\n")
cat("build_bias_correction.R finished at:", as.character(Sys.time()), "\n")
cat("Saved:", out_rds, "\n")
cat("Saved:", out_csv, "\n")
cat("Saved:", diag_csv, "\n")
cat("=============================================\n")
message("\n=== DONE: Bias correction factors ready ===")
