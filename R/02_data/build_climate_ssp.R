## =========================================================
## Project : r_project_tez
## Script  : R/02_data/build_climate_ssp.R
## Purpose : Build monthly M_climate(t) from CMIP6 tas for active SSP
##           M_climate scales importation pressure by global warming signal.
##
## Method  : M(y) = 1 + k * (dT(y) - dT(baseline)),  k = 0.13/°C
##           dT = area-weighted annual mean temperature anomaly
##           Baseline year = 2024  →  M(2024) = 1.0
##
## Model   : CNRM-CM6-1-HR (all SSP scenarios)
## Input   : DIR_CLIMATE_RAW/{tas}/*.nc  (set by SSP_SCENARIO in paths.R)
## Output  : DIR_PROCESSED_SSP/import_climate_multiplier_monthly.rds
##
## Reference: Cheng et al. (2023) eBioMedicine meta-analysis (54 studies):
##            RR = 1.13 per +1°C (95% CI: 1.11–1.16) → k ≈ 0.13
##            Sensitivity: k ∈ {0.10, 0.13, 0.20}
## =========================================================

source("R/01_setup/init.R")

## ---------------------------------------------------------
## 1) Parameters
## ---------------------------------------------------------
# k can be overridden via environment variable for sensitivity analysis
k_elasticity  <- as.numeric(Sys.getenv("K_ELASTICITY", unset = "0.13"))
                           # Default: Cheng et al. 2023 meta-analysis: RR=1.13/°C
                           # (54 studies, >4M dengue cases; eBioMedicine 98:104872)
                           # Sensitivity: K_ELASTICITY ∈ {0.10, 0.13, 0.20}
baseline_year <- 2024L
start_year    <- 2024L
end_year      <- 2100L

cat("\n=============================================\n")
cat("build_climate_ssp.R started at:", as.character(Sys.time()), "\n")
cat("SSP_SCENARIO  :", SSP_SCENARIO, "\n")
cat("DIR_CLIMATE_RAW:", DIR_CLIMATE_RAW, "\n")
cat("k_elasticity  :", k_elasticity, "\n")
cat("baseline_year :", baseline_year, "\n")
cat("=============================================\n\n")

## ---------------------------------------------------------
## 2) Locate tas NetCDF (SSP-specific)
## ---------------------------------------------------------
tas_dir <- file.path(DIR_CLIMATE_RAW, "tas")
if (!dir.exists(tas_dir)) {
  stop("Missing tas directory: ", tas_dir, call. = FALSE)
}

nc_files <- list.files(tas_dir, pattern = "\\.nc$", full.names = TRUE)
if (length(nc_files) == 0) {
  stop("No NetCDF files found in: ", tas_dir, call. = FALSE)
}

nc_path <- nc_files[1]
cat("Using NetCDF:", basename(nc_path), "\n\n")

## ---------------------------------------------------------
## 3) Read tas and compute area-weighted monthly global mean
## ---------------------------------------------------------
r  <- terra::rast(nc_path)
tt <- terra::time(r)
yy <- lubridate::year(tt)
mm <- lubridate::month(tt)

# Limit to analysis window
idx <- which(yy >= start_year & yy <= end_year)
r   <- r[[idx]]
tt  <- tt[idx]
yy  <- lubridate::year(tt)
mm  <- lubridate::month(tt)

# Area weights: cos(lat) for proper spherical averaging
lat <- terra::yFromRow(r, 1:nrow(r))
w   <- cos(pi * lat / 180)
W   <- terra::rast(r[[1]])
terra::values(W) <- rep(w, each = ncol(r))

# Monthly area-weighted mean temperature
gmean <- sapply(seq_len(terra::nlyr(r)), function(i) {
  xi <- terra::values(r[[i]])
  wi <- terra::values(W)
  ok <- is.finite(xi) & is.finite(wi)
  sum(xi[ok] * wi[ok]) / sum(wi[ok])
})

dfT <- tibble::tibble(
  date  = as.Date(tt),
  year  = yy,
  month = mm,
  Tglob = as.numeric(gmean)
) %>%
  dplyr::mutate(
    # Kelvin -> Celsius (CMIP6 tas is in K)
    T_C = dplyr::if_else(Tglob > 100, Tglob - 273.15, Tglob)
  )

cat("Temperature range:", round(min(dfT$T_C), 1), "–",
    round(max(dfT$T_C), 1), "°C\n")

## ---------------------------------------------------------
## 4) dT relative to baseline year mean
## ---------------------------------------------------------
T_baseline <- dfT %>%
  dplyr::filter(year == baseline_year) %>%
  dplyr::summarise(TC_base = mean(T_C, na.rm = TRUE)) %>%
  dplyr::pull(TC_base)

cat("T_baseline (", baseline_year, "):", round(T_baseline, 2), "°C\n\n")

dfM <- dfT %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(TC_year = mean(T_C, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    deltaT_centered = TC_year - T_baseline,
    M_climate_year  = 1 + k_elasticity * deltaT_centered
  )

## ---------------------------------------------------------
## 5) Expand to monthly (flat within year)
## ---------------------------------------------------------
dfM_monthly <- dfT %>%
  dplyr::select(year, month) %>%
  dplyr::distinct() %>%
  dplyr::left_join(
    dfM %>% dplyr::select(year, M_climate_year, deltaT_centered),
    by = "year"
  ) %>%
  dplyr::rename(M_climate = M_climate_year) %>%
  dplyr::arrange(year, month)

## ---------------------------------------------------------
## 6) Save to SSP-specific directory
## ---------------------------------------------------------
out_path <- file.path(DIR_PROCESSED_SSP, "import_climate_multiplier_monthly.rds")
saveRDS(dfM_monthly, out_path)

# Also save a CSV for inspection
out_csv <- file.path(DIR_PROCESSED_SSP, "import_climate_multiplier_monthly.csv")
readr::write_csv(dfM_monthly, out_csv)

cat("=============================================\n")
cat("build_climate_ssp.R finished\n")
cat("SSP scenario :", SSP_SCENARIO, "\n")
cat("Saved RDS    :", out_path, "\n")
cat("Saved CSV    :", out_csv, "\n")
cat("Years        :", min(dfM_monthly$year), "-", max(dfM_monthly$year), "\n")
cat("M(2024)      :", round(mean(dfM_monthly$M_climate[dfM_monthly$year == baseline_year]), 4), "\n")
cat("M(2050)      :", round(mean(dfM_monthly$M_climate[dfM_monthly$year == 2050]), 4), "\n")
cat("M(2075)      :", round(mean(dfM_monthly$M_climate[dfM_monthly$year == 2075]), 4), "\n")
cat("M(2100)      :", round(mean(dfM_monthly$M_climate[dfM_monthly$year == 2100]), 4), "\n")
cat("=============================================\n")

message("\n=== DONE (", SSP_SCENARIO, ") ===")
