# ============================================================
# OAT — TAM TİTİZLİK: 5 m-senaryosu × 3 SSP, n_rep=1000
# Gece boyu koşar (~10-12 saat). Ana koşumla aynı seed/n_mc → baz = ana tablo.
# Yalnızca m-senaryoları (beta/ip atlanır: OAT m-tablosu için gereksiz).
# ============================================================
m_grid <- c(0.5, 0.8, 1.0, 1.2, 2.0)

for (s in c("ssp126","ssp245","ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source("R/01_setup/init.R")
  
  # run_ctmc_spark'ı ANA koşumu tetiklemeden yükle (izole ortam)
  .E <- new.env(parent = globalenv())
  sys.source("R/03_models/ctmc_spark_monte_carlo.R", envir = .E)
  rcs <- .E$run_ctmc_spark
  
  out_root <- file.path(DIR_OUTPUT_SSP, "sensitivity", "ctmc_mc_rep1000")
  
  for (mv in m_grid) {
    scen <- sprintf("m_%03d", round(mv * 100))     # m_050, m_080, m_100, m_120, m_200
    cat(sprintf("\n>>> %s | %s (m=%.2f) | %s\n", s, scen, mv, format(Sys.time(), "%H:%M")))
    rcs(
      year_min = 2025, year_max = 2075, tau = 30L,
      infectious_period_days = 5, beta_vh = 0.3, beta_hv = 0.33,
      m = mv,                                       # ← değişen tek parametre
      include_import_in_q = FALSE, use_stochastic_EIP = TRUE,
      n_mc_eip = 2000L, n_rep = 1000L, seed = 123L,  # ← ana koşumla AYNI
      major_rule = list(type = "establishment"),
      out_dir = file.path(out_root, scen)
    )
  }
  cat(sprintf(">>> %s TAMAM\n", s))
}
cat("\n===== TÜM OAT KOŞUMLARI BİTTİ =====\n")

# ---------- oku + baz (m=1.0) ana tabloya eşit mi ----------
library(readr); library(dplyr); library(purrr)
did <- c(TUR.40.25_1="Kartal", TUR.10.4_1="Hopa", TUR.59.4_1="Fethiye",
         TUR.81.6_1="Zonguldak", TUR.39.3_1="Eğirdir")
sc  <- c(`m=0.5`="m_050", `m=0.8`="m_080", `m=1.0`="m_100", `m=1.2`="m_120", `m=2.0`="m_200")

for (s in c("ssp126","ssp245","ssp585")) {
  bd <- file.path("outputs", s, "sensitivity", "ctmc_mc_rep1000")
  t <- imap(sc, function(f, ml) {
    ff <- list.files(file.path(bd, f), "ctmc_spark_horizon.*rep1000.*csv$", full.names = TRUE)[1]
    read_csv(ff, show_col_types = FALSE) %>% transmute(district_id, !!ml := p_ge1_major_mean)
  }) %>% reduce(full_join, by = "district_id") %>%
    mutate(ilce = did[district_id], .after = district_id)
  
  mf <- list.files(file.path("outputs", s, "simulation"),
                   "ctmc_spark_horizon.*rep1000.*csv$", full.names = TRUE)
  main <- read_csv(mf[which.max(file.mtime(mf))], show_col_types = FALSE)
  t$ana_tablo <- main$p_ge1_major_mean[match(t$district_id, main$district_id)]
  
  cat("\n=========", s, "OAT n_rep=1000 + ana tablo =========\n")
  print(as.data.frame(t[, c("ilce","m=0.5","m=0.8","m=1.0","m=1.2","m=2.0","ana_tablo")]),
        digits = 4, row.names = FALSE)
  cat("KONTROL: 'm=1.0' ≈ 'ana_tablo' olmalı (n_rep=1000 → birebir).\n")
}