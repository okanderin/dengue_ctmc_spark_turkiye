# ============================================================
# TAM DUYARLILIK: 9 senaryo (m×5 + β±20% + ip±20%) × 3 SSP, n_rep=1000
# ~18-22 saat. Ana koşumla aynı seed/n_mc → baz = ana tablo.
# NOT: PRCC (Şekil 6.9.13) bundan ÜRETİLMEZ; o R0-tabanlı, ithalattan bağımsız, zaten doğru.
#      Beta/ip yalnızca TORNADO içindir.
# ============================================================
for (s in c("ssp126","ssp245","ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source("R/01_setup/init.R")
  
  # run_ctmc_spark'ı ANA koşumu tetiklemeden yükle
  .E <- new.env(parent = globalenv())
  sys.source("R/03_models/ctmc_spark_monte_carlo.R", envir = .E)
  assign("run_ctmc_spark", .E$run_ctmc_spark, envir = globalenv())
  
  # run_sensitivity_mc'yi yükle (dosya sonundaki otomatik n_rep=1000 koşumunu ATLA)
  .sl <- readLines("R/03_models/sensitivity_ctmc_mc.R")
  .ep <- grep("^# Entrypoint", .sl); .ep <- if (length(.ep)) .ep[1] else grep("run_sensitivity_mc\\(", .sl)[1]
  eval(parse(text = paste(.sl[seq_len(.ep - 1L)], collapse = "\n")), envir = globalenv())
  
  cat(sprintf("\n########## %s BAŞLADI %s ##########\n", s, format(Sys.time(), "%H:%M")))
  run_sensitivity_mc(n_rep = 1000L, n_mc_eip = 2000L, seed = 123L)   # 9 senaryo
  cat(sprintf("########## %s TAMAM %s ##########\n", s, format(Sys.time(), "%H:%M")))
}
cat("\n===== TÜM SENARYOLAR BİTTİ =====\n")

# ---------- oku: OAT (m) + tornado (β/ip) + baz kontrolü ----------
library(readr); library(dplyr); library(purrr)
did <- c(TUR.40.25_1="Kartal", TUR.10.4_1="Hopa", TUR.59.4_1="Fethiye",
         TUR.81.6_1="Zonguldak", TUR.39.3_1="Eğirdir")

for (s in c("ssp126","ssp245","ssp585")) {
  bd <- file.path("outputs", s, "sensitivity", "ctmc_mc")
  # mevcut tüm senaryo klasörlerini oku (base, m_*, beta_*, ip_*)
  folders <- list.dirs(bd, recursive = FALSE, full.names = FALSE)
  allv <- imap(setNames(folders, folders), function(f, nm) {
    ff <- list.files(file.path(bd, f), "ctmc_spark_horizon.*rep1000.*csv$", full.names = TRUE)
    if (!length(ff)) return(NULL)
    read_csv(ff[which.max(file.mtime(ff))], show_col_types = FALSE) %>%
      transmute(district_id, !!nm := p_ge1_major_mean)
  }) %>% compact() %>% reduce(full_join, by = "district_id") %>%
    mutate(ilce = did[district_id], .after = district_id)
  
  # ana tablo
  mf <- list.files(file.path("outputs", s, "simulation"),
                   "ctmc_spark_horizon.*rep1000.*csv$", full.names = TRUE)
  main <- read_csv(mf[which.max(file.mtime(mf))], show_col_types = FALSE)
  allv$ana_tablo <- main$p_ge1_major_mean[match(allv$district_id, main$district_id)]
  
  cat("\n=========", s, "— TÜM SENARYOLAR (n_rep=1000) =========\n")
  print(as.data.frame(allv[, setdiff(names(allv), "district_id")]), digits = 4, row.names = FALSE)
  cat("KONTROL: 'base' sütunu ≈ 'ana_tablo' olmalı (birebir).\n")
}