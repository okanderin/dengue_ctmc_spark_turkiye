# ============================================================
# (a) OAT'yi ZON ithalatıyla yeniden koş + (b) Regresyonu tek kaynaktan üret
# Çalışma dizini .Rproj kökü olmalı. Dış döngü atıl olduğundan OAT n_rep=1 ile
# koşulur (m=1,00 ana tabloya birebir eşit çıkar). Dakikalar sürer, saatler değil.
# ============================================================

## ---------- (a) OAT — üç SSP, n_rep=1 ----------
for (s in c("ssp126","ssp245","ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source("R/01_setup/init.R")                       # DIR_OUTPUT_SSP'yi bu SSP için türet
  
  # run_ctmc_spark'ı ANA koşumu TETİKLEMEDEN yükle (izole ortam → entrypoint pasif)
  .E <- new.env(parent = globalenv())
  sys.source("R/03_models/ctmc_spark_monte_carlo.R", envir = .E)
  assign("run_ctmc_spark", .E$run_ctmc_spark, envir = globalenv())
  
  # run_sensitivity_mc'yi yükle, dosya sonundaki (n_rep=1000) otomatik koşumu ATLA
  .sl <- readLines("R/03_models/sensitivity_ctmc_mc.R")
  .ep <- grep("^# Entrypoint", .sl); .ep <- if (length(.ep)) .ep[1] else grep("run_sensitivity_mc\\(", .sl)[1]
  eval(parse(text = paste(.sl[seq_len(.ep - 1L)], collapse = "\n")), envir = globalenv())
  
  run_sensitivity_mc(n_rep = 1L, n_mc_eip = 2000L, seed = 123L)   # atıl döngü → n_rep=1 yeterli
  cat(">>> OAT tamam:", s, "\n")
}

## ---------- OAT özetini oku + m=1,00 ana tabloya eşit mi? ----------
library(readr); library(dplyr)
for (s in c("ssp126","ssp245","ssp585")) {
  ss <- read_csv(file.path("outputs", s, "sensitivity", "ctmc_mc", "sensitivity_summary.csv"),
                 show_col_types = FALSE)
  cols <- intersect(c("district_id","m","beta_vh","ip_days","p_ge1_major_mean"), names(ss))
  cat("\n===", s, "— OAT özet ===\n"); print(as.data.frame(ss[, cols]), digits = 3)
  
  h <- read_csv(file.path("outputs", s, "simulation",
                          list.files(file.path("outputs", s, "simulation"),
                                     "ctmc_spark_horizon.*rep1000.*csv$", full.names = FALSE)[1]),
                show_col_types = FALSE)
  base <- ss %>% filter(abs(m - 1.0) < 1e-9, abs(beta_vh - 0.30) < 1e-9, ip_days == 5)
  cmp <- data.frame(
    ilce = c("Kartal","Hopa"),
    ana  = c(h$p_ge1_major_mean[h$district_id=="TUR.40.25_1"],
             h$p_ge1_major_mean[h$district_id=="TUR.10.4_1"]),
    OAT_m1 = c(base$p_ge1_major_mean[base$district_id=="TUR.40.25_1"],
               base$p_ge1_major_mean[base$district_id=="TUR.10.4_1"]))
  cat("m=1,00 vs ana tablo (eşit olmalı):\n"); print(cmp, digits = 4)
}

## ---------- (b) Regresyon — tek kaynak (6.10) ----------
Sys.setenv(SSP_SCENARIO = "ssp245"); source("R/01_setup/init.R")
source("R/04_results/6.10_pest_regression.R")     # mekanistik + Tobit + LMM → outputs/tables/
cat("\n--- Regresyon özeti ---\n")
print(read_csv("outputs/tables/topt_comparison.csv", show_col_types = FALSE))
print(read_csv("outputs/tables/cross_ssp_topt.csv",  show_col_types = FALSE))