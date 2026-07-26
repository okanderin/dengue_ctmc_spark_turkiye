# ==========================================================
# R/04_results/build_eip_sdlog_table.R
# EIP LogNormal parametrelerinin (mean EIP, sdlog, %5–%95) ilçe–ay tablosu
#
# - Kalibrasyon TAKLİT EDİLMEZ: parameter_functions.R içindeki
#   get_eip_logn_params() doğrudan çağrılır (tek doğruluk kaynağı).
# - Sıcaklık: climate_sentinel_monthly_wide_2015_2100.rds (tas, °C)
# - Tür eşleşmesi: sentinel_species.csv
# - Referans yıl: 2050 (SSP2-4.5), üreme sezonu Mayıs–Ekim
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
})

# ---- 0) Yollar (KENDİ PROJENİZE GÖRE AYARLAYIN) -----------
# Öneri: paths.R içindeki değişkenleri kullanın (ör. PROC_DIR, MODEL_DIR)
here_scenario <- "ssp245"
REF_YEAR      <- 2050L
MONTHS        <- 5:10   # Mayıs–Ekim

path_param_fns <- "R/03_models/parameter_functions.R"
path_traits    <- "data_processed/trait_params_aegypti.csv"      # varsayılan; tür bazlı aşağıda
path_traits_ae <- "data_processed/trait_params_aegypti.csv"
path_traits_al <- "data_processed/trait_params_albopictus.csv"
path_species   <- "data_processed/sentinel_species.csv"
path_climate   <- file.path("data_processed", here_scenario,
                            "climate_sentinel_monthly_wide_2015_2100.rds")

out_csv  <- file.path("outputs", here_scenario, "tables", "tbl_eip_sdlog_district_month.csv")
out_docx <- file.path("outputs", here_scenario, "tables", "tbl_eip_sdlog_district_month.docx")

# ---- 1) Fonksiyonları ve trait tablolarını yükle ----------
source(path_param_fns)  # get_eip_logn_params, eip_mean_of_T, read_trait_params ...

tab_ae <- read_trait_params(path_traits_ae)
tab_al <- read_trait_params(path_traits_al)
tab_by_species <- list(aegypti = tab_ae, albopictus = tab_al)

# ---- 2) İlçe adları ve tür eşleşmesi ----------------------
species_map <- read_csv(path_species, show_col_types = FALSE)

district_labels <- tibble::tribble(
  ~district_id,   ~ilce,       ~il,
  "TUR.59.4_1",   "Kartal",    "İstanbul",
  "TUR.81.6_1",   "Zonguldak", "Zonguldak",
  "TUR.10.4_1",   "Hopa",      "Artvin",
  "TUR.39.3_1",   "Fethiye",   "Muğla",
  "TUR.40.25_1",  "Eğirdir",   "Isparta"
)

month_tr <- c("1"="Ocak","2"="Şubat","3"="Mart","4"="Nisan","5"="Mayıs","6"="Haziran",
              "7"="Temmuz","8"="Ağustos","9"="Eylül","10"="Ekim","11"="Kasım","12"="Aralık")

# ---- 3) Sıcaklık verisini süz -----------------------------
clim <- readRDS(path_climate) %>%
  filter(year == REF_YEAR, month %in% MONTHS) %>%
  transmute(district_id, month = as.integer(month), tas = as.numeric(tas))

# ---- 4) Her ilçe–ay için EIP LogNormal parametrelerini üret -
# ÖNEMLİ: sdlog kalibrasyonu tür bazlı olduğundan, her tür için
# get_eip_logn_params çağrılmadan önce ilgili trait tablosu verilir.
# (parameter_functions.R oturum içinde kalibrasyonu cache'ler; tür
#  değişince force_recalibrate = TRUE ile yeniden kalibre ederiz.)

compute_row <- function(district_id, month, tas) {
  sp  <- species_map$species[species_map$district_id == district_id]
  tab <- tab_by_species[[sp]]

  # Tür değiştiğinde kalibrasyonu tazele (cache tek tür için geçerli)
  pars <- get_eip_logn_params(tas, tab, force_recalibrate = TRUE)

  if (!is.finite(pars$meanlog)) {
    return(tibble(mean_eip = Inf, sdlog = pars$sdlog,
                  p5 = Inf, p50 = Inf, p95 = Inf, species = sp))
  }
  q <- function(p) exp(pars$meanlog + pars$sdlog * qnorm(p))
  tibble(
    mean_eip = eip_mean_of_T(tas, tab),
    sdlog    = pars$sdlog,
    p5       = q(0.05),
    p50      = q(0.50),
    p95      = q(0.95),
    species  = sp
  )
}

tbl <- clim %>%
  mutate(res = pmap(list(district_id, month, tas), compute_row)) %>%
  unnest(res) %>%
  left_join(district_labels, by = "district_id") %>%
  mutate(
    Ay = month_tr[as.character(month)],
    across(c(tas, mean_eip, p5, p50, p95), ~round(.x, 1)),
    sdlog = round(sdlog, 3)
  ) %>%
  mutate(ilce = factor(ilce, levels = c("Kartal","Zonguldak","Hopa","Fethiye","Eğirdir"))) %>%
  arrange(ilce, month) %>%
  transmute(
    `İlçe` = ilce, `İl` = il, `Tür` = species, Ay,
    `T (°C)` = tas, `EIP ort. (gün)` = mean_eip, `sdlog(T)` = sdlog,
    `%5 (gün)` = p5, `Medyan` = p50, `%95 (gün)` = p95
  )

print(tbl, n = Inf)

# ---- 5) Kaydet -------------------------------------------
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write_csv(tbl, out_csv)

# İsteğe bağlı: Word tablosu (flextable)
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  ft <- flextable::flextable(tbl) |>
    flextable::merge_v(j = c("İlçe","İl","Tür")) |>
    flextable::theme_vanilla() |>
    flextable::fontsize(size = 9, part = "all") |>
    flextable::autofit()
  officer::read_docx() |>
    flextable::body_add_flextable(ft) |>
    print(target = out_docx)
  message("✓ DOCX yazıldı: ", out_docx)
}

message("✓ Tablo tamam: ", out_csv)
