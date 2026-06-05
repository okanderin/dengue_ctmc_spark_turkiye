# ==========================================================
# R/02_data/build_trait_params_species_files.R
#
# Build species-specific trait parameter files:
#   - trait_params_aegypti.csv
#   - trait_params_albopictus.csv
#
# Duzeltmeler (bu versiyonda):
#   1) mu_temp -> lifespan_temp  (lifespan parametreleri, mu_v = 1/lf)
#   2) lifespan_temp c degerleri duzeltildi (Mordecai 2017, lf olcegi)
#   3) vpd_ref ve k_vpd dolduruldu (Schmidt et al. 2018)
# ==========================================================

source("R/01_setup/init.R")

base_path <- file.path(DIR_PROCESSED, "trait_params.csv")
out_aeg   <- file.path(DIR_PROCESSED, "trait_params_aegypti.csv")
out_alb   <- file.path(DIR_PROCESSED, "trait_params_albopictus.csv")

if (!file.exists(base_path)) stop("Base trait file not found: ", base_path, call. = FALSE)

tab <- readr::read_csv(base_path, show_col_types = FALSE)

# ---- Schema asserts
req  <- c("trait","form","c","T0","Tm","unit","notes","source")
miss <- setdiff(req, names(tab))
if (length(miss) > 0) {
  stop("trait_params.csv missing columns: ", paste(miss, collapse=", "), call. = FALSE)
}

# ---- Required traits
need_traits <- c("a_biting", "eip_dev_rate", "lifespan_temp", "vpd_ref", "k_vpd")
miss_traits <- setdiff(need_traits, tab$trait)
if (length(miss_traits) > 0) {
  stop(
    "trait_params.csv missing trait rows: ", paste(miss_traits, collapse=", "),
    "\nHint: Run build_parameter_trait.R first to generate the template.",
    call. = FALSE
  )
}

# ==========================================================
# VPD parametreleri hakkinda not (her iki turde ayni):
#
#   vpd_ref = 1.0 kPa:
#     Referans VPD degeri; bu noktada nem carpani = 1 (etkisiz).
#     Schmidt et al. 2018 meta-analizinde Ae. aegypti icin en dusuk
#     mortalite ~27.5 degC / orta nem kosuluna denk gelen ~1 kPa
#     civarinda. Bu degerde exp(k*(vpd - vpd_ref)) = 1.
#
#   k_vpd = 0.5 per kPa:
#     VPD'nin mortaliteye duyarlilik katsayisi.
#     Brady et al. 2013 ve Schmidt et al. 2018 cercevesiyle uyumlu;
#     VPD 1 kPa arttikca mortalite yaklasik %65 artar (exp(0.5) ~ 1.65).
#     Duyarlilik analizi onerilir: k_vpd = 0.3, 0.5, 0.8.
# ==========================================================

VPD_REF <- 1.0   # kPa  (Schmidt et al. 2018)
K_VPD   <- 0.5   # /kPa (Brady et al. 2013; Schmidt et al. 2018)

# ==========================================================
# 1) Aedes aegypti
#
# Kaynak: Mordecai et al. 2017, PLOS NTD, Table S2
#   a_biting:      c=2.71e-4, T0=14.67, Tm=41.00  [1/gun]
#   eip_dev_rate:  c=1.04e-4, T0=11.50, Tm=38.97  [1/gun]
#   lifespan_temp: c=1.48e-1, T0=9.16,  Tm=37.73  [gun] <- lf, mu=1/lf
# ==========================================================

tab_aeg <- tab %>%
  dplyr::mutate(
    c = dplyr::case_when(
      trait == "a_biting"      ~ 2.71e-4,
      trait == "eip_dev_rate"  ~ 1.04e-4,
      trait == "lifespan_temp" ~ 1.48e-1,   # lf parametresi (gun olcegi)
      trait == "vpd_ref"       ~ VPD_REF,
      trait == "k_vpd"         ~ K_VPD,
      TRUE ~ c
    ),
    T0 = dplyr::case_when(
      trait == "a_biting"      ~ 14.67,
      trait == "eip_dev_rate"  ~ 11.50,
      trait == "lifespan_temp" ~ 9.16,
      TRUE ~ T0
    ),
    Tm = dplyr::case_when(
      trait == "a_biting"      ~ 41.00,
      trait == "eip_dev_rate"  ~ 38.97,
      trait == "lifespan_temp" ~ 37.73,
      TRUE ~ Tm
    ),
    form = dplyr::case_when(
      trait %in% c("a_biting","eip_dev_rate") ~ "briere",
      trait == "lifespan_temp"                ~ "quadratic",
      trait %in% c("vpd_ref","k_vpd")         ~ "scalar",
      TRUE ~ form
    ),
    unit = dplyr::case_when(
      trait %in% c("a_biting","eip_dev_rate") ~ "1/day",
      trait == "lifespan_temp"                ~ "days",
      trait == "vpd_ref"                      ~ "kPa",
      trait == "k_vpd"                        ~ "per_kPa",
      TRUE ~ unit
    ),
    notes = dplyr::case_when(
      trait == "lifespan_temp" ~
        "Adult lifespan lf(T) days; mu_v=1/lf in mu_v_of_TRH(); Ae. aegypti override",
      trait == "vpd_ref" ~
        "Reference VPD (kPa); Schmidt et al. 2018 / Brady et al. 2013",
      trait == "k_vpd" ~
        "VPD sensitivity of mortality (/kPa); Schmidt et al. 2018 / Brady et al. 2013",
      trait %in% c("a_biting","eip_dev_rate") ~
        paste0(notes, ifelse(is.na(notes) | notes=="", "", " | "),
               "Species override: Ae. aegypti"),
      TRUE ~ notes
    ),
    source = dplyr::case_when(
      trait %in% c("a_biting","eip_dev_rate","lifespan_temp") ~
        "Mordecai et al. 2017, PLOS NTD (Ae. aegypti thermal response)",
      trait %in% c("vpd_ref","k_vpd") ~
        "Schmidt et al. 2018 Parasites & Vectors; Brady et al. 2013 Parasites & Vectors",
      TRUE ~ source
    )
  )

# Sanity checks
check_traits_aeg <- c("a_biting","eip_dev_rate","lifespan_temp","vpd_ref","k_vpd")
if (any(is.na(tab_aeg$c[tab_aeg$trait %in% check_traits_aeg])))
  stop("Ae. aegypti: NA in c after override.", call.=FALSE)

# ==========================================================
# 2) Aedes albopictus
#
# Kaynak: Mordecai et al. 2017, PLOS NTD, Table S2
#   a_biting:      c=1.93e-4, T0=10.25, Tm=38.32  [1/gun]
#   eip_dev_rate:  c=1.09e-4, T0=10.39, Tm=43.05  [1/gun]
#   lifespan_temp: c=1.43e-1, T0=6.24,  Tm=38.25  [gun] <- lf, mu=1/lf
# ==========================================================

tab_alb <- tab %>%
  dplyr::mutate(
    c = dplyr::case_when(
      trait == "a_biting"      ~ 1.93e-4,
      trait == "eip_dev_rate"  ~ 1.09e-4,
      trait == "lifespan_temp" ~ 1.43e-1,   # lf parametresi (gun olcegi)
      trait == "vpd_ref"       ~ VPD_REF,
      trait == "k_vpd"         ~ K_VPD,
      TRUE ~ c
    ),
    T0 = dplyr::case_when(
      trait == "a_biting"      ~ 10.25,
      trait == "eip_dev_rate"  ~ 10.39,
      trait == "lifespan_temp" ~ 6.24,
      TRUE ~ T0
    ),
    Tm = dplyr::case_when(
      trait == "a_biting"      ~ 38.32,
      trait == "eip_dev_rate"  ~ 43.05,
      trait == "lifespan_temp" ~ 38.25,
      TRUE ~ Tm
    ),
    form = dplyr::case_when(
      trait %in% c("a_biting","eip_dev_rate") ~ "briere",
      trait == "lifespan_temp"                ~ "quadratic",
      trait %in% c("vpd_ref","k_vpd")         ~ "scalar",
      TRUE ~ form
    ),
    unit = dplyr::case_when(
      trait %in% c("a_biting","eip_dev_rate") ~ "1/day",
      trait == "lifespan_temp"                ~ "days",
      trait == "vpd_ref"                      ~ "kPa",
      trait == "k_vpd"                        ~ "per_kPa",
      TRUE ~ unit
    ),
    notes = dplyr::case_when(
      trait == "lifespan_temp" ~
        "Adult lifespan lf(T) days; mu_v=1/lf in mu_v_of_TRH(); Ae. albopictus override",
      trait == "vpd_ref" ~
        "Reference VPD (kPa); Schmidt et al. 2018 / Brady et al. 2013",
      trait == "k_vpd" ~
        "VPD sensitivity of mortality (/kPa); Schmidt et al. 2018 / Brady et al. 2013",
      trait %in% c("a_biting","eip_dev_rate") ~
        paste0(notes, ifelse(is.na(notes) | notes=="", "", " | "),
               "Species override: Ae. albopictus"),
      TRUE ~ notes
    ),
    source = dplyr::case_when(
      trait %in% c("a_biting","eip_dev_rate","lifespan_temp") ~
        "Mordecai et al. 2017, PLOS NTD (Ae. albopictus thermal response)",
      trait %in% c("vpd_ref","k_vpd") ~
        "Schmidt et al. 2018 Parasites & Vectors; Brady et al. 2013 Parasites & Vectors",
      TRUE ~ source
    )
  )

# Sanity checks
check_traits_alb <- c("a_biting","eip_dev_rate","lifespan_temp","vpd_ref","k_vpd")
if (any(is.na(tab_alb$c[tab_alb$trait %in% check_traits_alb])))
  stop("Ae. albopictus: NA in c after override.", call.=FALSE)

# ---- Write outputs
readr::write_csv(tab_aeg, out_aeg, na = "")
readr::write_csv(tab_alb, out_alb, na = "")

message("✅ Species-specific trait files written:")
message(" - ", out_aeg)
message(" - ", out_alb)
message("")
message("ℹ Parametreler:")
message("  Ae. aegypti    lifespan_temp: c=1.48e-1, T0=9.16,  Tm=37.73 [mu_v=1/lf]")
message("  Ae. albopictus lifespan_temp: c=1.43e-1, T0=6.24,  Tm=38.25 [mu_v=1/lf]")
message("  vpd_ref = ", VPD_REF, " kPa  |  k_vpd = ", K_VPD, " /kPa  (her iki tur)")