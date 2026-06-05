# ==========================================================
# build_travel_weights_static.R
# Purpose : Build STATIC (time-invariant) district travel weights
#           using TÜİK 2022 district-level accommodation arrivals
#           (ilce_konaklama_2022_yigm.xlsx)
#
# Key change (v2): district_factor = gelis_yabanci (foreign arrivals)
#   from TÜİK "Tesis Geliş Sayısı - Yabancı" column.
#   Previous version used pop_2024 * doluluk_pct which was a poor
#   proxy and overestimated Zonguldak/Kartal by 100-500x.
#
# Inputs:
#   - data_raw/population/turkstat/sentinel_baseline_pop.xlsx
#   - data_raw/population/turizm_verileri/ilce_konaklama_2022_yigm.xlsx
#
# Output:
#   - data_processed/travel_weights_static.csv
#   - data_processed/travel_weights_static.rds
#
# Columns in output:
#   district_id      : GADM code (e.g. TUR.40.25_1)
#   province_name    : province name
#   district_name    : district name
#   pop_2024         : 2024 population
#   gelis_yabanci    : TÜİK 2022 foreign arrivals (raw)
#   gelis_toplam     : TÜİK 2022 total arrivals (raw, for reference)
#   travel_weight    : normalized weight (sums to 1 across 5 districts)
# ==========================================================

# ---- 0) Bootstrap: detect ROOT + init -----------------------------------
suppressPackageStartupMessages({
  if (!requireNamespace("rprojroot", quietly = TRUE)) install.packages("rprojroot")
})

if (!exists("ROOT")) {
  ROOT <- tryCatch(
    rprojroot::find_root(rprojroot::has_file("r_project_tez.Rproj")),
    error = function(e) getwd()
  )
}

init_path <- file.path(ROOT, "R", "01_setup", "init.R")
if (!file.exists(init_path)) stop("init.R not found at: ", init_path, call. = FALSE)
source(init_path)  # defines DIR_POP_TURK, DIR_POP, DIR_PROCESSED, DIR_LOGS

# ---- 1) Packages ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(tidyr)
  if (!requireNamespace("janitor",  quietly = TRUE)) install.packages("janitor")
  if (!requireNamespace("stringi",  quietly = TRUE)) install.packages("stringi")
  library(janitor)
  library(stringi)
})

# ---- 2) Logging ----------------------------------------------------------
log_file <- file.path(DIR_LOGS, "build_travel_weights_static_log.txt")
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_line("Started build_travel_weights_static.R (v2 - TÜİK gelis_yabanci)")
log_line("ROOT = ", ROOT)

# ---- 3) Helpers ----------------------------------------------------------
std_text <- function(x) {
  x %>% as.character() %>% str_squish() %>% str_replace_all("\\s+", " ")
}

ascii_norm <- function(x) {
  x %>% str_to_lower() %>%
    str_replace_all("ş", "s") %>% str_replace_all("ı", "i") %>%
    str_replace_all("ğ", "g") %>% str_replace_all("ü", "u") %>%
    str_replace_all("ö", "o") %>% str_replace_all("ç", "c") %>%
    str_replace_all("\\s+", " ") %>% str_squish()
}

make_district_id <- function(province_name, district_name) {
  paste0(province_name, "_", district_name) %>%
    as.character() %>%
    stringr::str_squish() %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

# ---- 4) Baseline sentinel pop -------------------------------------------
baseline_path <- file.path(DIR_POP_TURK, "sentinel_baseline_pop.xlsx")
if (!file.exists(baseline_path)) stop("Missing: ", baseline_path, call. = FALSE)

baseline_raw <- readxl::read_excel(baseline_path) %>% janitor::clean_names()

col_prov <- intersect(names(baseline_raw), c("province", "il", "province_name"))[1]
col_dist <- intersect(names(baseline_raw), c("district", "ilce", "district_name"))[1]
col_pop  <- intersect(names(baseline_raw), c("pop_2024", "population", "pop", "nufus"))[1]

if (is.na(col_prov) || is.na(col_dist) || is.na(col_pop)) {
  stop("Baseline column mapping failed. Found: ",
       paste(names(baseline_raw), collapse = ", "), call. = FALSE)
}

baseline <- baseline_raw %>%
  transmute(
    province_name    = std_text(.data[[col_prov]]),
    district_name    = std_text(.data[[col_dist]]),
    pop_2024         = as.numeric(.data[[col_pop]])
  ) %>%
  mutate(
    district_id_slug = make_district_id(province_name, district_name),
    prov_key         = ascii_norm(province_name),
    dist_key         = ascii_norm(district_name)
  )

log_line("Baseline loaded: ", nrow(baseline), " districts")
log_line("  Districts: ", paste(baseline$district_name, collapse = ", "))

# ---- 5) TÜİK ilce_konaklama_2022 ----------------------------------------
# File layout (col indices, 0-based):
#   0: İL (province, Excel merged cells -> NAs after first row of group)
#   1: İLÇE (district)
#   2: Tesis Geliş YABANCI  <- district_factor for V_TR
#   3: Tesis Geliş YERLİ
#   4: Tesis Geliş TOPLAM
#   5-7:  Geceleme (YABANCI, YERLİ, TOPLAM)
#   8-10: Ortalama Kalış Süresi
#   11-13: Doluluk Oranı (%)
# First 3 rows are title + header rows, data starts row 4.

DIR_TOURISM    <- file.path(DIR_POP, "turizm_verileri")
konaklama_path <- file.path(DIR_TOURISM, "ilce_konaklama_2022_yigm.xlsx")

if (!file.exists(konaklama_path)) {
  stop("Missing TÜİK accommodation file: ", konaklama_path,
       "\nExpected: data_raw/population/turizm_verileri/ilce_konaklama_2022_yigm.xlsx",
       call. = FALSE)
}

kon_raw <- readxl::read_excel(konaklama_path, col_names = FALSE)
log_line("TÜİK raw rows: ", nrow(kon_raw))

# Drop 3 header rows, keep first 5 columns only
kon_data <- kon_raw %>%
  slice(-(1:3)) %>%
  dplyr::select(1:5)

colnames(kon_data) <- c("il_raw", "ilce_raw",
                         "gelis_yabanci", "gelis_yerli", "gelis_toplam")

# Forward-fill province (merged cells leave NAs)
kon_data <- kon_data %>%
  mutate(il_raw = dplyr::na_if(as.character(il_raw), "NA")) %>%
  tidyr::fill(il_raw, .direction = "down") %>%
  mutate(
    province_name = std_text(il_raw),
    district_name = std_text(as.character(ilce_raw)),
    gelis_yabanci = suppressWarnings(as.numeric(gelis_yabanci)),
    gelis_yerli   = suppressWarnings(as.numeric(gelis_yerli)),
    gelis_toplam  = suppressWarnings(as.numeric(gelis_toplam))
  ) %>%
  filter(
    !is.na(district_name), district_name != "NA", nchar(district_name) > 0,
    !is.na(gelis_yabanci)
  ) %>%
  mutate(
    prov_key = ascii_norm(province_name),
    dist_key = ascii_norm(district_name)
  ) %>%
  dplyr::select(prov_key, dist_key, gelis_yabanci, gelis_yerli, gelis_toplam)

log_line("TÜİK rows after cleaning: ", nrow(kon_data))

# ---- 6) Join baseline + TÜİK arrivals -----------------------------------
# Primary: ASCII province + district
joined <- baseline %>%
  dplyr::left_join(
    kon_data,
    by = c("prov_key", "dist_key")
  )

n_matched <- sum(!is.na(joined$gelis_yabanci))
log_line("Matched (province+district): ", n_matched, " / ", nrow(joined))

# Fallback: district name only (handles Zonguldak Merkez where il may differ)
if (any(is.na(joined$gelis_yabanci))) {
  unmatched <- joined %>% filter(is.na(gelis_yabanci)) %>%
    pull(district_name)
  log_line("Unmatched, trying district-only fallback: ",
           paste(unmatched, collapse = ", "))

  kon_dist <- kon_data %>%
    dplyr::select(dist_key, gelis_yabanci, gelis_yerli, gelis_toplam)

  joined <- joined %>%
    dplyr::left_join(kon_dist, by = "dist_key", suffix = c("", "_fb")) %>%
    dplyr::mutate(
      gelis_yabanci = dplyr::coalesce(gelis_yabanci, gelis_yabanci_fb),
      gelis_yerli   = dplyr::coalesce(gelis_yerli,   gelis_yerli_fb),
      gelis_toplam  = dplyr::coalesce(gelis_toplam,  gelis_toplam_fb)
    ) %>%
    dplyr::select(-dplyr::ends_with("_fb"))

  log_line("Matched after fallback: ", sum(!is.na(joined$gelis_yabanci)),
           " / ", nrow(joined))
}

# Hard stop if still unmatched
if (any(is.na(joined$gelis_yabanci))) {
  missing <- joined %>%
    filter(is.na(gelis_yabanci)) %>%
    dplyr::select(province_name, district_name, prov_key, dist_key)
  cat("Unmatched districts:\n"); print(missing)
  stop("Could not match all districts to TÜİK arrivals data.", call. = FALSE)
}

# ---- 7) Compute travel_weight -------------------------------------------
# V_TR_i = gelis_yabanci_i / sum(gelis_yabanci)
# Rationale: foreign accommodation arrivals is the direct proxy for
# import pressure. No power smoothing — raw proportions are used.

static_weights <- joined %>%
  mutate(
    travel_weight = gelis_yabanci / sum(gelis_yabanci, na.rm = TRUE)
  ) %>%
  dplyr::select(district_id_slug, province_name, district_name,
                pop_2024, gelis_yabanci, gelis_toplam, travel_weight)

# ---- 8) QA diagnostics --------------------------------------------------
log_line("--- Travel weight diagnostics (TÜİK 2022 gelis_yabanci) ---")
for (i in seq_len(nrow(static_weights))) {
  log_line(sprintf("  %-20s  gelis_yabanci=%7d  V_TR=%.4f",
                   static_weights$district_name[i],
                   as.integer(static_weights$gelis_yabanci[i]),
                   static_weights$travel_weight[i]))
}

s <- sum(static_weights$travel_weight, na.rm = TRUE)
if (!is.finite(s) || abs(s - 1) > 1e-10) {
  stop("QA failed: sum(travel_weight) = ", s, call. = FALSE)
}
log_line("sum(travel_weight) = ", format(s, scientific = FALSE), "  [OK]")

# ---- 9) GADM crosswalk --------------------------------------------------
gadm_crosswalk <- tibble::tribble(
  ~district_id_slug,   ~district_id,
  "artvin_hopa",        "TUR.10.4_1",
  "isparta_egirdir",    "TUR.39.3_1",
  "istanbul_kartal",    "TUR.40.25_1",
  "mugla_fethiye",      "TUR.59.4_1",
  "zonguldak_merkez",   "TUR.81.6_1"
)

static_weights <- static_weights %>%
  dplyr::left_join(gadm_crosswalk, by = "district_id_slug") %>%
  dplyr::mutate(
    district_id = dplyr::coalesce(district_id, district_id_slug)
  ) %>%
  dplyr::select(district_id, province_name, district_name,
                pop_2024, gelis_yabanci, gelis_toplam, travel_weight)

n_gadm <- sum(grepl("^TUR\\.", static_weights$district_id))
log_line("GADM id assigned: ", n_gadm, " / ", nrow(static_weights), " districts")

# ---- 10) Save -----------------------------------------------------------
out_csv <- file.path(DIR_PROCESSED, "travel_weights_static.csv")
out_rds <- file.path(DIR_PROCESSED, "travel_weights_static.rds")

readr::write_csv(static_weights, out_csv)
saveRDS(static_weights, out_rds)

log_line("Saved CSV: ", out_csv)
log_line("Saved RDS: ", out_rds)
log_line("Finished build_travel_weights_static.R (v2) successfully.")
