## =========================================================
## R/04_results/mc_validation_vs_analytic.R
##
## Purpose: Validate CTMC Monte Carlo implementation against
##          analytic birth-death (gambler's ruin) solution.
##
## Method:
##   1. Load ctmc_spark.R outputs (analytic P_est per district-month)
##   2. Load ctmc_spark_monte_carlo.R outputs (MC-averaged P_est)
##   3. Compare: |P_est_MC - P_est_analytic| < threshold
##
## Thesis reference:
##   "Model doğrulaması kapsamında, 5 ilçe için 50 yıllık
##    projeksiyon döneminde her yılın 12 ayı ayrı değerlendirilmiş
##    ve toplamda 3000 ilçe–yıl–ay kombinasyonu analiz edilmiştir."
##
## Outputs:
##   - mc_validation_results.csv (per district-month comparison)
##   - mc_validation_summary.csv (aggregate statistics)
##   - fig_mc_vs_analytic.png   (scatter plot)
##
## References:
##   - Allen (2017) Infect Dis Model — stochastic epidemic primer
##   - Keeling & Rohani (2008) — branching process theory
##   - Li (2018) — mathematical modeling validation
## =========================================================

source("R/01_setup/init.R")

SSP_SCENARIO <- Sys.getenv("SSP_SCENARIO", unset = "ssp245")

cat("\n=============================================\n")
cat("MC Validation vs Analytic — CTMC Spark Model\n")
cat("SSP:", SSP_SCENARIO, "\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("=============================================\n\n")

## =========================================================
## 1) Configuration
## =========================================================
THRESHOLD_ABS <- 0.02   # |P_mc - P_analytic| < 2%
THRESHOLD_REL <- 0.10   # |P_mc - P_analytic| / P_analytic < 10% (for P > 0.01)

## =========================================================
## 2) Load analytic outputs (ctmc_spark.R — single run)
## =========================================================
# Analytic run uses deterministic EIP → single P_est per district-month
# This is the "gold standard" gambler's ruin solution
analytic_dir <- file.path(DIR_OUTPUT_SSP, "model_results")
analytic_files <- list.files(analytic_dir,
                             pattern = "ctmc_spark_monthly.*\\.rds$",
                             full.names = TRUE)

if (length(analytic_files) == 0) {
  stop("No analytic (ctmc_spark.R) monthly output found in: ", analytic_dir,
       "\nRun ctmc_spark.R first with use_stochastic_EIP = FALSE.",
       call. = FALSE)
}

# Use the most recent file
analytic_path <- sort(analytic_files, decreasing = TRUE)[1]
cat("Analytic source:", basename(analytic_path), "\n")
df_analytic <- readRDS(analytic_path)

# Ensure required columns
req_cols_a <- c("district_id", "year", "month", "p_establishment",
                "lambda_local_i1", "temp_c")
missing_a <- setdiff(req_cols_a, names(df_analytic))
if (length(missing_a) > 0) {
  stop("Analytic output missing columns: ",
       paste(missing_a, collapse = ", "), call. = FALSE)
}

df_analytic <- df_analytic %>%
  dplyr::select(district_id, year, month,
                p_est_analytic = p_establishment,
                lambda_local_analytic = lambda_local_i1,
                temp_c) %>%
  dplyr::distinct()

cat("  Analytic rows:", nrow(df_analytic), "\n")
cat("  Districts:", dplyr::n_distinct(df_analytic$district_id), "\n")
cat("  Year range:", range(df_analytic$year), "\n\n")

## =========================================================
## 3) Load MC outputs (ctmc_spark_monte_carlo.R)
## =========================================================
mc_dir <- file.path(DIR_OUTPUT_SSP, "simulation")
mc_files <- list.files(mc_dir,
                       pattern = "ctmc_spark_monthly.*rep.*\\.rds$",
                       full.names = TRUE)

if (length(mc_files) == 0) {
  stop("No MC (ctmc_spark_monte_carlo.R) monthly output found in: ", mc_dir,
       "\nRun ctmc_spark_monte_carlo.R first.",
       call. = FALSE)
}

mc_path <- sort(mc_files, decreasing = TRUE)[1]
cat("MC source:", basename(mc_path), "\n")
df_mc <- readRDS(mc_path)

# MC outputs have _mean, _p2_5, _p97_5 columns
req_cols_m <- c("district_id", "year", "month",
                "p_establishment_mean", "p_establishment_p2_5",
                "p_establishment_p97_5")
missing_m <- setdiff(req_cols_m, names(df_mc))
if (length(missing_m) > 0) {
  stop("MC output missing columns: ",
       paste(missing_m, collapse = ", "), call. = FALSE)
}

df_mc <- df_mc %>%
  dplyr::select(district_id, year, month,
                p_est_mc_mean = p_establishment_mean,
                p_est_mc_p2_5 = p_establishment_p2_5,
                p_est_mc_p97_5 = p_establishment_p97_5,
                lambda_local_mc = lambda_local_i1_mean) %>%
  dplyr::distinct()

cat("  MC rows:", nrow(df_mc), "\n\n")

## =========================================================
## 4) Join and compute differences
## =========================================================
key <- c("district_id", "year", "month")

df_val <- df_analytic %>%
  dplyr::inner_join(df_mc, by = key)

cat("Matched rows:", nrow(df_val), "\n")
cat("Total combinations: 5 districts × 50 years × 12 months =",
    5 * 50 * 12, "(expected)\n\n")

df_val <- df_val %>%
  dplyr::mutate(
    # Absolute difference (Denklem 17)
    abs_diff = abs(p_est_mc_mean - p_est_analytic),

    # Relative difference (for non-trivial P_est)
    rel_diff = dplyr::if_else(
      p_est_analytic > 0.001,
      abs_diff / p_est_analytic,
      NA_real_
    ),

    # Within threshold?
    within_abs_threshold = abs_diff < THRESHOLD_ABS,
    within_rel_threshold = is.na(rel_diff) | rel_diff < THRESHOLD_REL,

    # Is analytic value within MC confidence interval?
    analytic_in_mc_ci = (p_est_analytic >= p_est_mc_p2_5) &
                        (p_est_analytic <= p_est_mc_p97_5),

    # Theoretical MC standard error (Denklem 18)
    # SE = sqrt(p(1-p)/n) where n = n_rep (1000)
    se_theoretical = sqrt(p_est_mc_mean * (1 - p_est_mc_mean) / 1000)
  )

## =========================================================
## 5) Summary statistics
## =========================================================
cat("--- Validation Summary ---\n\n")

# Overall
n_total <- nrow(df_val)
n_pass_abs <- sum(df_val$within_abs_threshold)
n_pass_rel <- sum(df_val$within_rel_threshold, na.rm = TRUE)
n_in_ci <- sum(df_val$analytic_in_mc_ci)

cat("Total comparisons:", n_total, "\n")
cat("Within absolute threshold (<", THRESHOLD_ABS, "):",
    n_pass_abs, "/", n_total,
    "(", round(n_pass_abs / n_total * 100, 1), "%)\n")
cat("Analytic within MC 95% CI:",
    n_in_ci, "/", n_total,
    "(", round(n_in_ci / n_total * 100, 1), "%)\n")

cat("\nAbsolute difference statistics:\n")
cat("  Mean:   ", signif(mean(df_val$abs_diff), 4), "\n")
cat("  Median: ", signif(median(df_val$abs_diff), 4), "\n")
cat("  Max:    ", signif(max(df_val$abs_diff), 4), "\n")
cat("  P95:    ", signif(quantile(df_val$abs_diff, 0.95), 4), "\n")

# By district
district_summary <- df_val %>%
  dplyr::group_by(district_id) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_abs_diff = mean(abs_diff),
    max_abs_diff = max(abs_diff),
    pct_within_threshold = mean(within_abs_threshold) * 100,
    pct_in_ci = mean(analytic_in_mc_ci) * 100,

    mean_p_est_analytic = mean(p_est_analytic),
    mean_p_est_mc = mean(p_est_mc_mean),

    # Pearson correlation
    cor_r = cor(p_est_analytic, p_est_mc_mean, use = "complete.obs"),

    .groups = "drop"
  )

cat("\n--- Per-District Summary ---\n")
print(district_summary, n = Inf, width = 120)

## =========================================================
## 6) Note on expected differences (stochastic vs deterministic EIP)
## =========================================================
# IMPORTANT: ctmc_spark.R with use_stochastic_EIP = TRUE already
# computes E[P_est(EIP_i)] via inner MC. The outer MC in
# ctmc_spark_monte_carlo.R repeats this process n_rep times.
#
# Therefore, the comparison is:
#   Analytic: E[P_est(EIP)] from single inner MC (2000 draws)
#   MC:       mean over n_rep outer repetitions of E[P_est(EIP)]
#
# Differences arise from:
#   (a) Different random seeds → different EIP draws
#   (b) Finite MC variance
#
# If ctmc_spark.R runs with use_stochastic_EIP = FALSE (deterministic),
# then the comparison is truly analytic vs MC, and differences should
# be zero up to MC sampling noise.
#
# For thesis: report both comparisons if available.

## =========================================================
## 7) Save outputs
## =========================================================
out_dir <- file.path(DIR_OUTPUT_SSP, "diagnostics", "mc_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Detailed results
readr::write_csv(df_val,
                 file.path(out_dir, "mc_validation_results.csv"))

# Summary
readr::write_csv(district_summary,
                 file.path(out_dir, "mc_validation_summary.csv"))

## =========================================================
## 8) Diagnostic plot
## =========================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Scatter: analytic vs MC
  p1 <- ggplot(df_val, aes(x = p_est_analytic, y = p_est_mc_mean)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = district_id), alpha = 0.3, size = 0.8) +
    scale_x_log10(labels = scales::scientific) +
    scale_y_log10(labels = scales::scientific) +
    labs(
      title = "CTMC Spark: Analitik vs Monte Carlo P_est",
      subtitle = paste0("n = ", nrow(df_val),
                        " | Mutlak fark < ", THRESHOLD_ABS, ": ",
                        round(n_pass_abs / n_total * 100, 1), "%"),
      x = expression(P[est]^{analitik}),
      y = expression(P[est]^{MC}),
      color = "İlçe"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(out_dir, "fig_mc_vs_analytic.png"),
         p1, width = 8, height = 6, dpi = 150)

  # Histogram of absolute differences
  p2 <- ggplot(df_val, aes(x = abs_diff)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = THRESHOLD_ABS, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    annotate("text", x = THRESHOLD_ABS, y = Inf,
             label = paste0("Eşik = ", THRESHOLD_ABS),
             hjust = -0.1, vjust = 2, color = "red", size = 3.5) +
    labs(
      title = "P_est Mutlak Fark Dağılımı (Analitik − MC)",
      x = expression("|"~P[est]^{analitik} - P[est]^{MC}~"|"),
      y = "Frekans"
    ) +
    theme_minimal(base_size = 11)

  ggsave(file.path(out_dir, "fig_abs_diff_histogram.png"),
         p2, width = 7, height = 4.5, dpi = 150)

  cat("\nPlots saved to:", out_dir, "\n")
}


## =========================================================
## 9) Pass/Fail verdict
## =========================================================
pass_rate <- n_pass_abs / n_total

cat("\n=============================================\n")
if (pass_rate >= 0.95) {
  cat("VALIDATION PASSED\n")
  cat("  ", round(pass_rate * 100, 1),
      "% of comparisons within absolute threshold\n")
} else if (pass_rate >= 0.90) {
  cat("VALIDATION MARGINAL\n")
  cat("  ", round(pass_rate * 100, 1),
      "% of comparisons within threshold (target: 95%)\n")
  cat("  Review per-district summary for outliers.\n")
} else {
  cat("VALIDATION FAILED\n")
  cat("  Only ", round(pass_rate * 100, 1),
      "% of comparisons within threshold\n")
  cat("  Investigate discrepancies before proceeding.\n")
}
cat("=============================================\n")

cat("\nOutputs:", out_dir, "\n")
cat("Finished:", as.character(Sys.time()), "\n")
