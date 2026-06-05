# ==========================================================
# Build sentinel_species.csv
# Species assignment per sentinel district
# ==========================================================

source("R/01_setup/init.R")

# ---- Output path
out_path <- file.path(DIR_PROCESSED, "sentinel_species.csv")

# ---- Sentinel species definition (based on current evidence)
# NOTE:
# TUR.10.4_1  -> Artvin (Hopa sentinel)        -> Aedes aegypti
# TUR.81.6_1  -> Zonguldak (Merkez sentinel)   -> Aedes aegypti
# Others      -> Aedes albopictus

sentinel_species <- tibble::tibble(
  district_id = c(
    "TUR.10.4_1",  # Artvin (Hopa)
    "TUR.81.6_1",  # Zonguldak (Merkez)
    "TUR.39.3_1",  # Isparta (Eğirdir)
    "TUR.40.25_1", # Istanbul(Kartal)
    "TUR.59.4_1"   # Mugla (Fethiye)
  ),
  species = c(
    "aegypti",
    "aegypti",
    "albopictus",
    "albopictus",
    "albopictus"
  )
)

# ---- Basic validation
if (any(is.na(sentinel_species$district_id))) {
  stop("district_id contains NA.", call. = FALSE)
}
if (any(!sentinel_species$species %in% c("aegypti","albopictus"))) {
  stop("species must be one of: aegypti, albopictus", call. = FALSE)
}
if (any(duplicated(sentinel_species$district_id))) {
  stop("Duplicate district_id detected.", call. = FALSE)
}

# ---- Write CSV
readr::write_csv(sentinel_species, out_path)

message("✅ sentinel_species.csv written to:\n", out_path)