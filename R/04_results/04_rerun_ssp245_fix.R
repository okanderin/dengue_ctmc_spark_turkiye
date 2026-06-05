## =========================================================
## 04_rerun_ssp245_fix.R
## SSP2-4.5 Lambda_import ≈ 0 hatasını düzeltmek için
## tam yeniden çalıştırma scripti.
##
## KÖK NEDEN:
##   ctmc_spark_monte_carlo.R satır 566 hardcoded olarak
##   "importation_pressure_monthly_2024_2100.rds" dosyasını kullanıyordu.
##   Bu dosya SSP2-4.5'te eski (pre-GBD) versiyondan kalmıştı.
##   Doğru dosya: "importation_pressure_monthly_2025_2075.rds" (v3, GBD-weighted).
##
## ÖN KOŞUL:
##   Düzeltilmiş ctmc_spark_monte_carlo.R ve run_sensitivity_importation.R
##   dosyalarını ilgili dizinlere kopyalamış olmalısınız:
##     - R/03_models/ctmc_spark_monte_carlo.R
##     - R/02_data/run_sensitivity_importation.R
##
## TAHMİNİ SÜRE: SSP başına ~45-90 dk (n_rep=1000, n_mc_eip=2000)
## =========================================================

cat("\n", strrep("=", 60), "\n")
cat("04_rerun_ssp245_fix.R — SSP2-4.5 Lambda_import düzeltmesi\n")
cat("Başlangıç:", as.character(Sys.time()), "\n")
cat(strrep("=", 60), "\n\n")

## ---- ADIM 0: Ön kontrol ----
cat(">>> ADIM 0: Ön kontrol <<<\n")

# ctmc_spark_monte_carlo.R'deki düzeltmeyi doğrula
mc_path <- here::here("R", "03_models", "ctmc_spark_monte_carlo.R")
mc_lines <- readLines(mc_path)
has_fix <- any(grepl("importation_pressure_monthly_2025_2075\\.rds", mc_lines))

if (!has_fix) {
  stop(
    "HATA: ctmc_spark_monte_carlo.R henüz düzeltilmemiş!\n",
    "Düzeltilmiş dosyayı R/03_models/ dizinine kopyalayın.\n",
    "Beklenen: 'importation_pressure_monthly_2025_2075.rds' referansı",
    call. = FALSE
  )
}
cat("  ✓ ctmc_spark_monte_carlo.R düzeltmesi doğrulandı.\n")

# İthalat dosyalarını kontrol et
for (s in c("ssp126", "ssp245", "ssp585")) {
  imp_path <- file.path(here::here("data_processed", s),
                        "importation_pressure_monthly_2025_2075.rds")
  if (!file.exists(imp_path)) {
    stop("HATA: İthalat dosyası bulunamadı: ", imp_path,
         "\nbuild_importation_pressure_monthly.R'yi önce çalıştırın.", call. = FALSE)
  }
  imp <- readRDS(imp_path)
  cat(sprintf("  ✓ %s: lambda_import mean=%.4f, max=%.4f, rows=%d\n",
              s, mean(imp$lambda_import, na.rm = TRUE),
              max(imp$lambda_import, na.rm = TRUE), nrow(imp)))
}

## ---- ADIM 1: SSP2-4.5 simülasyonu (zorunlu) ----
cat("\n>>> ADIM 1: SSP2-4.5 simülasyonu yeniden çalıştırılıyor <<<\n")
cat("  Bu adım ~45-90 dakika sürebilir.\n")

Sys.setenv(SSP_SCENARIO = "ssp245")
source(here::here("R", "01_setup", "init.R"))
source(here::here("R", "03_models", "ctmc_spark_monte_carlo.R"))

t0 <- Sys.time()
run_ctmc_spark()
cat("  SSP2-4.5 tamamlandı:", format(Sys.time() - t0), "\n")

# Doğrulama
hor_files <- list.files(here::here("outputs", "ssp245", "simulation"),
                        "ctmc_spark_horizon.*rep.*\\.csv$", full.names = TRUE)
h <- read.csv(sort(hor_files, decreasing = TRUE)[1])
cat("\n  SSP2-4.5 Lambda_import doğrulaması:\n")
print(h[, c("district_id", "Lambda_import")])

if (max(h$Lambda_import) < 100) {
  warning("Lambda_import hâlâ düşük — düzeltme etkili olmamış olabilir!")
} else {
  cat("  ✓ Lambda_import değerleri makul aralıkta.\n")
}

## ---- ADIM 2: SSP1-2.6 ve SSP5-8.5 (isteğe bağlı, tutarlılık için) ----
cat("\n>>> ADIM 2: SSP1-2.6 ve SSP5-8.5 (isteğe bağlı) <<<\n")
cat("  Bu senaryolarda zaten doğru dosya kullanılmış olabilir.\n")
cat("  Yine de tutarlılık için yeniden çalıştırmak önerilir.\n")

rerun_others <- TRUE  # FALSE yaparak atlayabilirsiniz

if (rerun_others) {
  for (s in c("ssp126", "ssp585")) {
    cat("\n  --- ", toupper(s), " ---\n")
    Sys.setenv(SSP_SCENARIO = s)
    source(here::here("R", "01_setup", "init.R"))
    source(here::here("R", "03_models", "ctmc_spark_monte_carlo.R"))
    
    t1 <- Sys.time()
    run_ctmc_spark()
    cat("  ", s, "tamamlandı:", format(Sys.time() - t1), "\n")
  }
}

## ---- ADIM 3: İthalat duyarlılık analizi (SSP2-4.5) ----
cat("\n>>> ADIM 3: İthalat duyarlılık analizi <<<\n")

Sys.setenv(SSP_SCENARIO = "ssp245")
source(here::here("R", "01_setup", "init.R"))
tryCatch({
  source(here::here("R", "02_data", "run_sensitivity_importation.R"))
  cat("  ✓ İthalat duyarlılık analizi tamamlandı.\n")
}, error = function(e) {
  cat("  ⚠ İthalat duyarlılık hatası:", e$message, "\n")
})

## ---- ADIM 4: Grafik ve tablo üretimi (01_generate) ----
cat("\n>>> ADIM 4: Grafik ve tablo üretimi <<<\n")

for (s in c("ssp126", "ssp245", "ssp585")) {
  Sys.setenv(SSP_SCENARIO = s)
  source(here::here("R", "04_results", "01_generate_ssp_outputs.R"))
}

## ---- ADIM 5: bulgular.Rmd knit ----
cat("\n>>> ADIM 5: bulgular.Rmd knit ediliyor <<<\n")

rmarkdown::render(here::here("R", "04_results", "bulgular.Rmd"))

cat("\n", strrep("=", 60), "\n")
cat("TAMAMLANDI:", as.character(Sys.time()), "\n")
cat("Çıktılar: outputs/ssp*/simulation/ ve outputs/ssp*/figures/\n")
cat(strrep("=", 60), "\n")
