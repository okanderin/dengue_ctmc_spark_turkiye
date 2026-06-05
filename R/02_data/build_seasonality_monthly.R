## =========================================================
## Project : r_project_tez
## Script  : R/02_data/build_seasonality_monthly.R
## Purpose : Build monthly seasonality multiplier S_m from
##           TÜİK/EGM visitor departure statistics.
##
## VERSION 3 — TÜİK "Çıkış Yapan Ziyaretçi Sayısı" (2024-2025)
##
## Previous version (v2) derived S_m from Türkiye-wide hotel
## occupancy rates (YIGM tesis doluluk oranı). While occupancy
## correlates with travel volume, it is a proxy that is bounded
## by fixed bed capacity and does not distinguish foreign from
## domestic visitors.
##
## This version uses EGM passport-based monthly departure counts
## (TÜİK Turizm İstatistikleri, Tablo 24100), which directly
## measure the quantity of interest: the seasonal variation in
## international visitor flow through Türkiye.
##
## Source : TÜİK "Çıkış Yapan Ziyaretçi Sayısı, Aylık"
##          (departures_monthly_2024_2025.csv)
##          https://data.tuik.gov.tr — Table 24100
##
## Normalization : S_m = monthly_avg / mean(monthly_avg_1..12)
##                 => mean(S_m) = 1, sum(S_m) = 12
##
## Outputs : data_processed/seasonality_monthly.csv/.rds
## Logs    : logs/build_seasonality_monthly_log.txt
## =========================================================

source("R/01_setup/init.R")

log_file <- file.path(DIR_LOGS, "build_seasonality_monthly_log.txt")
sink(log_file, split = TRUE)
on.exit({ try(sink(), silent = TRUE) }, add = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

cat("=============================================\n")
cat("build_seasonality_monthly.R (v3 — TÜİK departures)\n")
cat("Started at :", as.character(Sys.time()), "\n")
cat("ROOT       :", ROOT, "\n")
cat("=============================================\n\n")

## =========================================================
## 1) Read TÜİK monthly departure data
## =========================================================
## The CSV has a non-standard pivot format from TÜİK MEDAS.
## We hard-code the parsed values for robustness and
## reproducibility. Raw file is retained for audit trail.
##
## Source: TÜİK "Çıkış Yapan Ziyaretçi Sayısı, Aylık"
##         (departures_monthly_2024_2025.csv)
## Definition: "Yabancı" = T.C. pasaportu taşımayan,
##   yurt dışında ikamet eden ziyaretçiler (turist + günübirlikçi)
##   + yurt dışında ikamet eden T.C. vatandaşları.
##   EGM pasaport kayıtlarına dayalı tam sayım verisidir.

dep_2024 <- c(
  3245893, 2684623, 3079574, 4040974, 5585634, 6471276,
  7395701, 8661295, 7149583, 6723819, 4018247, 3175828
)

dep_2025 <- c(
  3355818, 2782758, 2982576, 4242307, 5618724, 6551137,
  7403653, 8912005, 7324078, 7079267, 4255110, 3409623
)

# Use 2024-2025 average to smooth year-to-year variation
dep_avg <- (dep_2024 + dep_2025) / 2

cat("--- Monthly departure counts (TÜİK, 2024-2025 average) ---\n")
month_names_tr <- c("Ocak","Subat","Mart","Nisan","Mayis","Haziran",
                     "Temmuz","Agustos","Eylul","Ekim","Kasim","Aralik")
for (m in 1:12) {
  cat(sprintf(
    "  %02d %-8s : %s\n",
    m,
    month_names_tr[m],
    format(dep_avg[m], big.mark = ",", scientific = FALSE)
  ))
}
cat(sprintf(
  "  Annual total: %s\n",
  format(sum(dep_avg), big.mark = ",", scientific = FALSE)
))

cat(sprintf(
  "  Peak (Aug)  : %s\n",
  format(max(dep_avg), big.mark = ",", scientific = FALSE)
))

cat(sprintf(
  "  Trough (Feb): %s\n",
  format(min(dep_avg), big.mark = ",", scientific = FALSE)
))

## =========================================================
## 2) Normalize to mean = 1
## =========================================================
## S_m = monthly_avg / mean(monthly_avg)
## This ensures sum(S_m) = 12 and mean(S_m) = 1.
## Interpretation: S_m > 1 means above-average travel month;
##                 S_m < 1 means below-average.

seasonality <- tibble(
  month = 1L:12L,
  departure_count = dep_avg,
  seasonal_multiplier = dep_avg / mean(dep_avg)
)

cat("\n--- Seasonality multiplier S_m (mean=1) ---\n")
print(seasonality %>% dplyr::select(month, seasonal_multiplier), n = 12)
cat("Mean(S_m) =", mean(seasonality$seasonal_multiplier), "\n")

rng <- range(seasonality$seasonal_multiplier)
cat("Range(S_m) =", round(rng[1], 4), "-", round(rng[2], 4), "\n")
cat("Peak/Trough ratio =", round(rng[2] / rng[1], 2), "\n")

## =========================================================
## 3) Compare with previous occupancy-based S_m (diagnostic)
## =========================================================
sm_occ_v2 <- c(0.628, 0.655, 0.652, 0.801, 1.103, 1.249,
               1.434, 1.534, 1.321, 1.244, 0.725, 0.654)

cat("\n--- Comparison: departure-based vs occupancy-based S_m ---\n")
cat(sprintf("  %2s  %8s  %8s  %7s\n", "Mo", "S_m(dep)", "S_m(occ)", "Diff%"))
for (m in 1:12) {
  diff_pct <- (seasonality$seasonal_multiplier[m] - sm_occ_v2[m]) /
    sm_occ_v2[m] * 100
  cat(sprintf("  %2d  %8.4f  %8.4f  %+6.1f%%\n",
              m, seasonality$seasonal_multiplier[m], sm_occ_v2[m], diff_pct))
}
cat("  Occupancy-based peak/trough:", round(max(sm_occ_v2)/min(sm_occ_v2), 2), "\n")
cat("  Departure-based peak/trough:", round(rng[2]/rng[1], 2), "\n")
cat("  Note: occupancy proxy underestimates seasonality amplitude\n")
cat("        due to fixed bed-capacity ceiling effect.\n")

## =========================================================
## 4) Save outputs
## =========================================================
out <- seasonality %>%
  dplyr::select(month, seasonal_multiplier)

out_csv <- file.path(DIR_PROCESSED, "seasonality_monthly.csv")
out_rds <- file.path(DIR_PROCESSED, "seasonality_monthly.rds")

readr::write_csv(out, out_csv)
saveRDS(out, out_rds)

cat("\nSaved:", out_csv, "\n")
cat("Saved:", out_rds, "\n")
cat("\n=== DONE (v3 — TÜİK departures) ===\n")
