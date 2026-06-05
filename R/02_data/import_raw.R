## =========================================================
## Project : r_project_tez
## Script  : R/02_data/import_raw.R
## Purpose : Import CMIP6 NetCDFs for active SSP scenario
##           and aggregate to sentinel districts (ADM2)
## Model   : CNRM-CM6-1-HR (all scenarios)
## Inputs  : DIR_CLIMATE_RAW (NetCDF folders per variable)
##           DIR_SHP (district boundary shapefile)
## Output  : DIR_INTERIM_SSP/climate_district_monthly_2015_2100.rds
## Log     : DIR_LOGS/import_raw_log.txt
## =========================================================

## ---- 0) Initialize project environment (packages + paths + global options)
source("R/01_setup/init.R")

## ---------------------------------------------------------
## 1) Logging (console + file)
## ---------------------------------------------------------
log_file <- file.path(DIR_LOGS, paste0("import_raw_", SSP_SCENARIO, "_log.txt"))

sink(log_file, split = TRUE)
on.exit({ try(sink(), silent = TRUE) }, add = TRUE)

cat("=============================================\n")
cat("import_raw.R started at :", as.character(Sys.time()), "\n")
cat("Project ROOT            :", ROOT, "\n")
cat("SSP_SCENARIO            :", SSP_SCENARIO, "\n")
cat("DIR_CLIMATE_RAW         :", DIR_CLIMATE_RAW, "\n")
cat("DIR_SHP                 :", DIR_SHP, "\n")
cat("DIR_INTERIM_SSP         :", DIR_INTERIM_SSP, "\n")
cat("=============================================\n\n")

## ---------------------------------------------------------
## 2) Read GADM ADM2 shapefile explicitly
## ---------------------------------------------------------
shp_path <- file.path(DIR_SHP, "gadm41_TUR_2.shp")
if (!file.exists(shp_path)) {
  stop("ADM2 shapefile not found: ", shp_path, call. = FALSE)
}
message("Using ADM2 shapefile: ", shp_path)

districts_raw <- sf::st_read(shp_path, quiet = TRUE)

## ---- Standardize geometry column name to "geometry"
geom_col <- attr(districts_raw, "sf_column")
if (is.null(geom_col) || !geom_col %in% names(districts_raw)) {
  stop("Could not detect sf geometry column (sf_column attribute missing).", call. = FALSE)
}
if (geom_col != "geometry") {
  names(districts_raw)[names(districts_raw) == geom_col] <- "geometry"
  attr(districts_raw, "sf_column") <- "geometry"
}

## ---- Ensure expected columns exist
needed <- c("NAME_1", "NAME_2", "GID_2")
missing_cols <- setdiff(needed, names(districts_raw))
if (length(missing_cols) > 0) {
  stop(
    "Expected GADM columns missing: ", paste(missing_cols, collapse = ", "), "\n",
    "Available columns: ", paste(names(districts_raw), collapse = ", "),
    call. = FALSE
  )
}

## ---------------------------------------------------------
## 3) Sentinel table (Province + District) and matching
## ---------------------------------------------------------
sentinel <- tibble::tibble(
  prov = c("Artvin",   "Zinguldak", "Istanbul", "Isparta", "Mugla"),
  dist = c("Hopa",     "Merkez",    "Kartal",  "Eğirdir", "Fethiye")
)

# Robust normalization: case/space + Turkish chars -> ASCII (stable join)
norm_tr <- function(x) {
  x <- stringr::str_squish(x)
  x <- stringr::str_to_lower(x)
  stringi::stri_trans_general(x, "Latin-ASCII")
}

districts_raw <- districts_raw %>%
  dplyr::mutate(
    prov_norm = norm_tr(NAME_1),
    dist_norm = norm_tr(NAME_2)
  )

sentinel <- sentinel %>%
  dplyr::mutate(
    prov_norm = norm_tr(prov),
    dist_norm = norm_tr(dist)
  )

districts <- districts_raw %>%
  dplyr::inner_join(
    sentinel %>% dplyr::select(prov_norm, dist_norm),
    by = c("prov_norm", "dist_norm")
  ) %>%
  dplyr::mutate(
    province_name = NAME_1,
    district_name = NAME_2,
    district_id   = as.character(GID_2)
  ) %>%
  dplyr::select(district_id, province_name, district_name, geometry)

## ---- Sanity checks (must match all sentinel rows)
cat("\n--- Sentinel matching report ---\n")
cat("Requested sentinel units:\n")
print(sentinel %>% dplyr::select(prov, dist))

cat("\nMatched units in GADM ADM2:\n")
print(districts %>% sf::st_drop_geometry() %>%
        dplyr::select(province_name, district_name, district_id))

if (nrow(districts) != nrow(sentinel)) {
  have <- districts_raw %>% dplyr::distinct(prov_norm, dist_norm)
  missing <- dplyr::anti_join(sentinel, have, by = c("prov_norm", "dist_norm")) %>%
    dplyr::select(prov, dist)

  cat("\nERROR: Some sentinel units were not found in GADM ADM2:\n")
  print(missing)

  stop("Sentinel matching failed. Fix spellings in sentinel table.", call. = FALSE)
}

## ---- CRS standardization: EPSG:4326
districts <- sf::st_make_valid(districts)
districts <- sf::st_transform(districts, 4326)

message("\nDistrict polygons loaded (sentinel): ", nrow(districts))
message("District CRS                     : EPSG:4326")
message("District IDs                     : ",
        paste(districts$district_id, collapse = ", "))

## ---------------------------------------------------------
## 4) Helper functions
## ---------------------------------------------------------
list_nc_files <- function(var_dir) {
  if (!dir.exists(var_dir)) stop("Missing variable folder: ", var_dir, call. = FALSE)
  files <- list.files(var_dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(files) == 0) stop("No NetCDF found in: ", var_dir, call. = FALSE)
  sort(files)
}

read_var_raster <- function(nc_files) {
  terra::rast(nc_files)
}

convert_units <- function(r, varname) {
  # tas, tasmin, tasmax: Kelvin -> Celsius
  # pr: kg m-2 s-1 -> mm/day
  # hurs: assumed % (no conversion)
  if (varname %in% c("tas", "tasmin", "tasmax")) {
    r <- r - 273.15
  } else if (varname == "pr") {
    r <- r * 86400
  }
  r
}

get_layer_dates <- function(r) {
  tt <- terra::time(r)
  if (is.null(tt) || all(is.na(tt))) {
    stop(
      "NetCDF time could not be parsed by terra::time().\n",
      "Check CF-compliance (time units/calendar) in NetCDF metadata.",
      call. = FALSE
    )
  }
  as.Date(tt)
}

district_means <- function(r, districts_sf) {
  exactextractr::exact_extract(r, districts_sf, fun = "mean", progress = FALSE)
}

## ---------------------------------------------------------
## 5) Variables to import (from paths.R CLIMATE_VARS)
## ---------------------------------------------------------
vars <- CLIMATE_VARS   # named vector: canonical_name = folder_name
cat("\nVariables to import:", paste(names(vars), collapse = ", "), "\n")
cat("Folder names       :", paste(vars, collapse = ", "), "\n\n")

## ---------------------------------------------------------
## 6) Main loop
## ---------------------------------------------------------
all_long <- list()

for (i in seq_along(vars)) {
  canonical <- names(vars)[i]     # internal name (e.g. "hurs")
  folder    <- vars[i]            # folder name in data_raw (e.g. "hurs")

  message("\n>>> Importing variable: ", canonical, " (folder: ", folder, ")")

  v_dir    <- file.path(DIR_CLIMATE_RAW, folder)
  nc_files <- list_nc_files(v_dir)

  message("NetCDF files (n=", length(nc_files), "):")
  for (f in nc_files) message("  - ", basename(f))

  r <- read_var_raster(nc_files)

  ## Align raster CRS to EPSG:4326 (districts CRS)
  if (!grepl("EPSG:4326|WGS 84|longlat", terra::crs(r), ignore.case = TRUE)) {
    message("Raster CRS is not EPSG:4326; projecting to EPSG:4326 ...")
    r <- terra::project(r, "EPSG:4326")
  }

  ## Unit conversion
  r <- convert_units(r, canonical)

  ## Dates
  dates <- get_layer_dates(r)
  message("Layer count : ", terra::nlyr(r))
  message("Date range  : ", as.character(min(dates)),
          " -> ", as.character(max(dates)))

  ## District means (sentinel only)
  m <- district_means(r, districts)

  ## Long format
  m_long <- as.data.frame(m) %>%
    dplyr::mutate(district_id = districts$district_id) %>%
    tidyr::pivot_longer(
      cols = -district_id,
      names_to = "layer",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      layer_index = as.integer(stringr::str_extract(layer, "\\d+$"))
    )

  if (any(is.na(m_long$layer_index))) {
    bad <- unique(m_long$layer[is.na(m_long$layer_index)])
    stop(
      "Could not parse layer index from some layer names.\n",
      "Examples: ", paste(head(bad, 10), collapse = ", "), "\n",
      "Inspect raster layer names and adjust regex if needed.",
      call. = FALSE
    )
  }

  m_long <- m_long %>%
    dplyr::mutate(
      date = dates[layer_index],
      year = lubridate::year(date),
      month = lubridate::month(date),
      variable = canonical
    ) %>%
    dplyr::left_join(
      districts %>% sf::st_drop_geometry() %>%
        dplyr::select(district_id, province_name, district_name),
      by = "district_id"
    ) %>%
    dplyr::select(district_id, province_name, district_name,
                  year, month, date, variable, value) %>%
    dplyr::arrange(district_id, date)

  all_long[[canonical]] <- m_long
  message("Rows produced for ", canonical, " : ", nrow(m_long))
}

climate_district_monthly <- dplyr::bind_rows(all_long) %>%
  dplyr::filter(year >= 2015, year <= 2100) %>%
  dplyr::arrange(district_id, date, variable)

## ---------------------------------------------------------
## 7) Save interim output (SSP-specific directory)
## ---------------------------------------------------------
out_file <- file.path(DIR_INTERIM_SSP,
                      "climate_district_monthly_2015_2100.rds")
saveRDS(climate_district_monthly, out_file)

cat("\n=============================================\n")
cat("import_raw.R finished at :", as.character(Sys.time()), "\n")
cat("SSP scenario             :", SSP_SCENARIO, "\n")
cat("Boundary file            :", shp_path, "\n")
cat("Output file              :", out_file, "\n")
cat("Total rows               :", nrow(climate_district_monthly), "\n")
cat("Districts                :",
    dplyr::n_distinct(climate_district_monthly$district_id), "\n")
cat("Date range               :",
    as.character(min(climate_district_monthly$date)),
    "->",
    as.character(max(climate_district_monthly$date)), "\n")
cat("Variables                :",
    paste(unique(climate_district_monthly$variable), collapse = ", "),
    "\n")
cat("=============================================\n")

message("\n=== DONE (", SSP_SCENARIO, ") ===")
message("Saved to: ", out_file)
