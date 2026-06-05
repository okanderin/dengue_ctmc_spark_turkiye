# ==========================================================
# R/01_setup/paths.R
# Central path definitions — SSP-scenario aware
#
# SSP_SCENARIO is read from environment variable or defaults to "ssp245".
# Set before sourcing init.R:
#   Sys.setenv(SSP_SCENARIO = "ssp126")
#   source("R/01_setup/init.R")
#
# Climate model: CNRM-CM6-1-HR (all SSP scenarios)
# ==========================================================

library(here)

# --- Active SSP scenario ---
SSP_SCENARIO <- Sys.getenv("SSP_SCENARIO", unset = "ssp245")
stopifnot(SSP_SCENARIO %in% c("ssp126", "ssp245", "ssp585"))
message("\u2139 Active SSP scenario: ", SSP_SCENARIO)

# --- Root directories ---
DIR_RAW       <- here("data_raw")
DIR_INTERIM   <- here("data_interim")
DIR_PROCESSED <- here("data_processed")
DIR_OUTPUT    <- here("outputs")
DIR_LOGS      <- here("logs")

# --- SSP-specific directories ---
DIR_CLIMATE_RAW   <- here(DIR_RAW, "climate", "cmip6", SSP_SCENARIO)
DIR_INTERIM_SSP   <- here(DIR_INTERIM, SSP_SCENARIO)
DIR_PROCESSED_SSP <- here(DIR_PROCESSED, SSP_SCENARIO)
DIR_OUTPUT_SSP    <- here(DIR_OUTPUT, SSP_SCENARIO)

# --- Historical / observation data (shared, not SSP-specific) ---
DIR_CLIMATE_HIST  <- here(DIR_RAW, "climate", "cmip6_hist_1981_2014")
DIR_ERA5          <- here(DIR_RAW, "climate", "era_5_land_1981_2014")

# --- Geography ---
DIR_SHP <- here(DIR_RAW, "geography", "districts_shapefile")

# --- Population data ---
DIR_POP       <- here(DIR_RAW, "population")
DIR_POP_TURK  <- here(DIR_RAW, "population", "turkstat")
DIR_POP_PROJ  <- here(DIR_RAW, "population", "projections")

# --- SSP-independent processed data (shared across scenarios) ---
DIR_PROCESSED_SHARED <- DIR_PROCESSED

# --- Cross-scenario comparison outputs ---
DIR_OUTPUT_CROSS <- here(DIR_OUTPUT, "cross_scenario")

# --- Standard input file paths (scenario-dependent) ---
fp_climate_interim <- here(DIR_INTERIM_SSP,
                           "climate_district_monthly_2015_2100.rds")
fp_climate_long    <- here(DIR_PROCESSED_SSP,
                           "climate_sentinel_monthly_long_2015_2100.rds")
fp_climate_wide    <- here(DIR_PROCESSED_SSP,
                           "climate_sentinel_monthly_wide_2015_2100.rds")
fp_import          <- here(DIR_PROCESSED_SSP,
                           "importation_pressure_monthly_2025_2075.rds")
fp_seasonality     <- here(DIR_PROCESSED_SSP,
                           "seasonality_monthly.rds")

# --- Standard input file paths (scenario-independent) ---
fp_population  <- here(DIR_PROCESSED_SHARED,
                       "population_sentinel_yearly_2015_2100.rds")
fp_species     <- here(DIR_PROCESSED_SHARED,
                       "sentinel_species.csv")
fp_travel      <- here(DIR_PROCESSED_SHARED,
                       "travel_weights_static.rds")
fp_trait       <- here(DIR_PROCESSED_SHARED,
                       "trait_params.csv")

# --- Climate variable mapping ---
# CNRM-CM6-1-HR uses "hurs" for near-surface relative humidity.
# This map translates canonical variable names to folder names in data_raw.
CLIMATE_VARS <- c(
  tas    = "tas",
  tasmin = "tasmin",
  tasmax = "tasmax",
  pr     = "pr",
  hurs   = "hurs"
)

# Internal canonical name used by the model pipeline (downstream scripts
# refer to the variable as "hurs"; the old name "hur" is no longer valid).
VAR_HUMIDITY <- "hurs"

# --- Create dirs if needed ---
for (d in c(DIR_INTERIM_SSP, DIR_PROCESSED_SSP, DIR_OUTPUT_SSP, DIR_LOGS,
            here(DIR_OUTPUT_SSP, "diagnostics"),
            here(DIR_OUTPUT_SSP, "figures"),
            here(DIR_OUTPUT_SSP, "model_results"),
            here(DIR_OUTPUT_SSP, "sensitivity"),
            here(DIR_OUTPUT_SSP, "simulation"),
            here(DIR_OUTPUT_SSP, "tables"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
