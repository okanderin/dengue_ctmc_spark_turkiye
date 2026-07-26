## =========================================================
## R/04_results/00_results_setup.R
## Paylaşımlı sabitler, tema ve yardımcı fonksiyonlar
## Tüm 04_results scriptleri tarafından source() edilir
## =========================================================

# ---- Paketler ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(scales)
  library(patchwork)
  library(kableExtra)
  library(flextable)
  library(officer)
  library(ggrepel)
})

# ---- Sabitler ----
AY_TR <- c("Oca","Şub","Mar","Nis","May","Haz",
           "Tem","Ağu","Eyl","Eki","Kas","Ara")
names(AY_TR) <- month.abb

DISTRICT_LABELS <- c(
  "TUR.10.4_1"  = "Hopa (Artvin)",
  "TUR.39.3_1"  = "Eğirdir (Isparta)",
  "TUR.40.25_1" = "Kartal (İstanbul)",
  "TUR.59.4_1"  = "Fethiye (Muğla)",
  "TUR.81.6_1"  = "Zonguldak"
)

DISTRICT_ORDER <- c(
  "TUR.40.25_1", "TUR.59.4_1", "TUR.81.6_1",
  "TUR.10.4_1",  "TUR.39.3_1"
)

# Renk paleti — district_id anahtarlı (aes(colour = district_id) için)
COL_DISTRICT <- c(
  "TUR.40.25_1" = "#E63946",
  "TUR.59.4_1"  = "#457B9D",
  "TUR.81.6_1"  = "#2A9D8F",
  "TUR.10.4_1"  = "#E9C46A",
  "TUR.39.3_1"  = "#264653"
)

# Renk paleti — district_label anahtarlı (aes(colour = district_label) için)
COL_DISTRICT_LABEL <- setNames(COL_DISTRICT, DISTRICT_LABELS[names(COL_DISTRICT)])

SSP_LABELS <- c(
  ssp126 = "SSP1-2.6",
  ssp245 = "SSP2-4.5",
  ssp585 = "SSP5-8.5"
)

# Renk paleti — ssp kod anahtarlı (aes(colour = ssp) için)
COL_SSP <- c(
  ssp126 = "#2A9D8F",
  ssp245 = "#E9C46A",
  ssp585 = "#E63946"
)

# Renk paleti — ssp_label anahtarlı (aes(colour = ssp_label) için)
COL_SSP_LABEL <- setNames(COL_SSP, SSP_LABELS[names(COL_SSP)])

SCEN_TR <- c(
  base          = "Baz (m=1.00, β=0.30, IP=5g)",
  m_050         = "m = 0.50 (düşük yoğunluk)",
  m_080         = "m = 0.80 (orta yoğunluk)",
  m_120         = "m = 1.20 (yüksek yoğunluk)",
  m_200         = "m = 2.00 (çok yüksek yoğunluk)",
  beta_minus_20 = "β−%20 (β=0.24)",
  beta_plus_20  = "β+%20 (β=0.36)",
  ip_plus_20    = "IP+%20 (IP=6g, γ↓)",
  ip_minus_20   = "IP−%20 (IP=4g, γ↑)"
)

# ---- Tema ----
theme_thesis <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      legend.position  = "bottom",
      legend.key.size  = unit(0.4, "cm"),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey40"),
      plot.caption     = element_text(size = base_size - 2, colour = "grey50", hjust = 0)
    )
}

# ---- Yardımcı fonksiyonlar ----
load_ssp_data <- function(ssp, base_dir = here("outputs")) {
  ssp_dir <- file.path(base_dir, ssp)
  sim_dir <- file.path(ssp_dir, "simulation")

  # Find most recent MC output
  monthly_files <- list.files(sim_dir, "ctmc_spark_monthly.*rep.*\\.csv$",
                              full.names = TRUE)
  yearly_files  <- list.files(sim_dir, "ctmc_spark_yearly.*rep.*\\.csv$",
                              full.names = TRUE)
  horizon_files <- list.files(sim_dir, "ctmc_spark_horizon.*rep.*\\.csv$",
                              full.names = TRUE)

  if (length(monthly_files) == 0)
    stop("No MC monthly output for ", ssp, " in ", sim_dir)

  monthly <- read_csv(monthly_files[which.max(file.mtime(monthly_files))],
                      show_col_types = FALSE) %>%
    mutate(
      ssp = ssp,
      district_id    = factor(district_id, levels = DISTRICT_ORDER),
      district_label = DISTRICT_LABELS[as.character(district_id)],
      month_name     = factor(AY_TR[month], levels = AY_TR)
    )

  yearly <- read_csv(yearly_files[which.max(file.mtime(yearly_files))],
                     show_col_types = FALSE) %>%
    mutate(
      ssp = ssp,
      district_id    = factor(district_id, levels = DISTRICT_ORDER),
      district_label = DISTRICT_LABELS[as.character(district_id)]
    )

  horizon <- read_csv(horizon_files[which.max(file.mtime(horizon_files))],
                      show_col_types = FALSE) %>%
    mutate(
      ssp = ssp,
      district_id    = factor(district_id, levels = DISTRICT_ORDER),
      district_label = DISTRICT_LABELS[as.character(district_id)]
    )

  list(monthly = monthly, yearly = yearly, horizon = horizon)
}

# Load importation sensitivity
load_sens_imp <- function(ssp, base_dir = here("outputs")) {
  path <- file.path(base_dir, ssp, "sensitivity", "importation",
                    paste0("sensitivity_summary_", ssp, ".csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(ssp = ssp)
}

# Load model (CTMC) sensitivity
load_sens_mc <- function(ssp, base_dir = here("outputs")) {
  path <- file.path(base_dir, ssp, "sensitivity", "ctmc_mc",
                    "sensitivity_summary.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      ssp = ssp,
      district_id    = factor(district_id, levels = DISTRICT_ORDER),
      district_label = DISTRICT_LABELS[as.character(district_id)]
    )
}

# Load country contributions
load_country_contrib <- function(ssp, base_dir = here("data_processed")) {
  path <- file.path(base_dir, ssp,
                    "importation_country_contributions.csv")
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>% mutate(ssp = ssp)
}

message("✓ 00_results_setup.R loaded")
