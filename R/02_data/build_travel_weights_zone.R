# ==========================================================
# build_travel_weights_zone.R
# Purpose : Build ZONE (province-level) travel weights for the 5 sentinel
#           districts. Each sentinel is treated as the CLIMATE-ARCHETYPE
#           representative of its PROVINCE; the importation weight uses the
#           PROVINCE total foreign accommodation arrivals (Σ over all
#           districts of that province), not the sentinel district alone.
#
# Rationale: Kartal itself is residential (low foreign arrivals), but it
#   represents the İstanbul climate zone, whose other districts receive
#   heavy tourism. Aggregating to the province makes the "representative
#   receptor" interpretation CONSISTENT WITH THE WEIGHTS (unlike relabelling
#   the district-level weights, which would contradict the arithmetic).
#
# Drop-in: sets `gelis_yabanci` = province total, so the Method-1
#   importation script (V_in_day = gelis_yabanci × S_m/365) becomes
#   ZONE-level automatically when pointed at this file.
#
# λ note: π_weighted and S_m are district-INDEPENDENT (joined by year/month
#   only), so λ_import ∝ A across districts. Hence, per sentinel,
#     Λ_zone / Λ_district = A_province / A_district
#   (printed in the comparison block below — no need to re-run importation
#   to see how the ranking shifts).
#
# Inputs:
#   - data_raw/population/turkstat/sentinel_baseline_pop.xlsx
#   - data_raw/population/turizm_verileri/ilce_konaklama_2022_yigm.xlsx
# Output:
#   - data_processed/travel_weights_zone.csv / .rds
#
# Output columns:
#   district_id, province_name, district_name, pop_2024
#   gelis_yabanci_district   : sentinel district's own foreign arrivals (ref)
#   gelis_yabanci_province   : PROVINCE total foreign arrivals (zone A)
#   gelis_yabanci            : = province total (drop-in A_d for Method 1)
#   gelis_toplam_province    : province total arrivals (ref)
#   n_districts_in_province  : how many districts aggregated (diagnostic)
#   travel_weight_district   : district-based weight (sums to 1 over 5) [ref]
#   travel_weight            : ZONE weight (sums to 1 over 5) [active]
#   weight_basis             : "province"
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
log_file <- file.path(DIR_LOGS, "build_travel_weights_zone_log.txt")
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_line("Started build_travel_weights_zone.R (province-level A)")
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

log_line("Baseline loaded: ", nrow(baseline), " sentinels")
log_line("  Sentinels: ", paste(baseline$district_name, collapse = ", "))

# ---- 5) TÜİK ilce_konaklama_2022 (district accommodation, all provinces) -
# Layout (0-based): 0 İL | 1 İLÇE | 2 Geliş YABANCI | 3 YERLİ | 4 TOPLAM | ...
# First 3 rows are title/header; data starts row 4.
DIR_TOURISM    <- file.path(DIR_POP, "turizm_verileri")
konaklama_path <- file.path(DIR_TOURISM, "ilce_konaklama_2022_yigm.xlsx")

if (!file.exists(konaklama_path)) {
  stop("Missing TÜİK accommodation file: ", konaklama_path,
       "\nExpected: data_raw/population/turizm_verileri/ilce_konaklama_2022_yigm.xlsx",
       call. = FALSE)
}

kon_raw <- readxl::read_excel(konaklama_path, col_names = FALSE)
log_line("TÜİK raw rows: ", nrow(kon_raw))

kon_data <- kon_raw %>%
  slice(-(1:3)) %>%
  dplyr::select(1:5)

colnames(kon_data) <- c("il_raw", "ilce_raw",
                        "gelis_yabanci", "gelis_yerli", "gelis_toplam")

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
  dplyr::select(province_name, prov_key, dist_key,
                gelis_yabanci, gelis_yerli, gelis_toplam)

log_line("TÜİK rows after cleaning: ", nrow(kon_data))

# ---- 6) DISTRICT value per sentinel (for reference / comparison) --------
# Primary: province+district; fallback: district-only (as in v2).
dist_primary <- baseline %>%
  dplyr::left_join(
    kon_data %>% dplyr::select(prov_key, dist_key,
                              gelis_yabanci_district = gelis_yabanci),
    by = c("prov_key", "dist_key")
  )

if (any(is.na(dist_primary$gelis_yabanci_district))) {
  kon_dist_only <- kon_data %>%
    dplyr::group_by(dist_key) %>%
    dplyr::summarise(gelis_yabanci_fb = sum(gelis_yabanci, na.rm = TRUE),
                     .groups = "drop")
  dist_primary <- dist_primary %>%
    dplyr::left_join(kon_dist_only, by = "dist_key") %>%
    dplyr::mutate(
      gelis_yabanci_district = dplyr::coalesce(gelis_yabanci_district,
                                               gelis_yabanci_fb)
    ) %>%
    dplyr::select(-gelis_yabanci_fb)
}

if (any(is.na(dist_primary$gelis_yabanci_district))) {
  stop("Bazi sentinel ilcelerinin ilce-duzeyi gelisi eslesmedi.", call. = FALSE)
}

# ---- 7) PROVINCE (zone) aggregate: Σ foreign arrivals over the province -
kon_province <- kon_data %>%
  dplyr::group_by(prov_key) %>%
  dplyr::summarise(
    province_name           = dplyr::first(province_name),
    gelis_yabanci_province  = sum(gelis_yabanci, na.rm = TRUE),
    gelis_toplam_province   = sum(gelis_toplam,  na.rm = TRUE),
    n_districts_in_province = dplyr::n(),
    .groups = "drop"
  )

zone <- dist_primary %>%
  dplyr::left_join(
    kon_province %>% dplyr::select(prov_key,
                                  gelis_yabanci_province,
                                  gelis_toplam_province,
                                  n_districts_in_province),
    by = "prov_key"
  )

if (any(is.na(zone$gelis_yabanci_province))) {
  miss <- zone %>% filter(is.na(gelis_yabanci_province)) %>%
    dplyr::select(province_name, district_name, prov_key)
  cat("Il eslesmeyen sentineller:\n"); print(miss)
  stop("Bazi sentinel illeri TUIK il-toplamina eslesmedi.", call. = FALSE)
}

# ---- 8) Weights: district-based (ref) and ZONE-based (active) -----------
zone <- zone %>%
  mutate(
    travel_weight_district = gelis_yabanci_district /
      sum(gelis_yabanci_district, na.rm = TRUE),
    travel_weight          = gelis_yabanci_province /
      sum(gelis_yabanci_province, na.rm = TRUE),   # ACTIVE (zone)
    gelis_yabanci          = gelis_yabanci_province, # drop-in A_d for Method 1
    weight_basis           = "province"
  )

# ---- 9) COMPARISON: district vs zone (λ scales linearly with A) ---------
log_line("--- District vs ZONE comparison (Lambda_zone/Lambda_district = A_prov/A_dist) ---")
cmp <- zone %>%
  transmute(
    district_name, province_name,
    A_district = round(gelis_yabanci_district),
    A_province = round(gelis_yabanci_province),
    ratio_prov_dist = round(gelis_yabanci_province / gelis_yabanci_district, 1),
    w_district = round(travel_weight_district, 4),
    w_zone     = round(travel_weight, 4)
  ) %>%
  dplyr::arrange(dplyr::desc(w_zone))
print(as.data.frame(cmp), row.names = FALSE)

cat("\nIthalat agirlik siralamasi:\n")
cat("  ILCE tabani: ",
    paste(zone$district_name[order(-zone$travel_weight_district)], collapse = " > "), "\n")
cat("  ZON  tabani: ",
    paste(zone$district_name[order(-zone$travel_weight)], collapse = " > "), "\n")
cat("\nNOT: Lambda ilcelerarasi A ile dogrusal oldugundan, yukaridaki\n",
    "ratio_prov_dist sutunu her sentinelin Lambda'sinin ilce->zon\n",
    "yeniden-olcekleme faktorudur. Nihai P_ufuk siralamasi Lambda x P_est\n",
    "carpimindan gelir; P_est degismedigi icin recombine yeterli (MC yok).\n", sep = "")

# ---- 10) GADM crosswalk --------------------------------------------------
gadm_crosswalk <- tibble::tribble(
  ~district_id_slug,   ~district_id,
  "artvin_hopa",        "TUR.10.4_1",
  "isparta_egirdir",    "TUR.39.3_1",
  "istanbul_kartal",    "TUR.40.25_1",
  "mugla_fethiye",      "TUR.59.4_1",
  "zonguldak_merkez",   "TUR.81.6_1"
)

out <- zone %>%
  dplyr::left_join(gadm_crosswalk, by = "district_id_slug") %>%
  dplyr::mutate(district_id = dplyr::coalesce(district_id, district_id_slug)) %>%
  dplyr::select(district_id, province_name, district_name, pop_2024,
                gelis_yabanci_district, gelis_yabanci_province,
                gelis_yabanci, gelis_toplam_province,
                n_districts_in_province,
                travel_weight_district, travel_weight, weight_basis)

# ---- 11) QA + Save -------------------------------------------------------
s <- sum(out$travel_weight, na.rm = TRUE)
if (!is.finite(s) || abs(s - 1) > 1e-10) {
  stop("QA failed: sum(travel_weight zone) = ", s, call. = FALSE)
}
log_line("sum(travel_weight zone) = ", format(s, scientific = FALSE), "  [OK]")

out_csv <- file.path(DIR_PROCESSED, "travel_weights_zone.csv")
out_rds <- file.path(DIR_PROCESSED, "travel_weights_zone.rds")
readr::write_csv(out, out_csv)
saveRDS(out, out_rds)

log_line("Saved CSV: ", out_csv)
log_line("Saved RDS: ", out_rds)
log_line("Finished build_travel_weights_zone.R successfully.")

cat("\n>>> Zon agirliklari ile ithalat uretmek icin:\n")
cat("    build_importation_pressure_monthly.R (v4) icindeki tw_path'i\n")
cat("    'travel_weights_zone.rds' yapin (ya da dosyayi kopyalayin;\n")
cat("    orijinal travel_weights_static.rds'i once yedekleyin).\n")
