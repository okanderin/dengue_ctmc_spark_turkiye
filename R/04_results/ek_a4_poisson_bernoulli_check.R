# =============================================================================
# ek_a4_poisson_bernoulli_check.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Reproducibly regenerate Appendix Table A.4 ("P_ufuk analizinin Poisson
#   Binomial Dagilimi (Tumleyen Carpimi) ile Yakinsamasi").
#
#   The table compares, for each district x SSP scenario, the EXACT
#   heterogeneous-Bernoulli complement product (Equation 12) against the
#   Poisson cumulative-hazard approximation (Equation A.14), and reports the
#   second-order Taylor error bound Sum(p^2)/2.
#
#   Until now this table was produced ad hoc (never scripted), which caused the
#   Section 6.5.1 prose to drift out of sync with the table (stale 0.484/0.482
#   vs. canonical 0.494/0.4925). This script makes the table a deterministic
#   function of the monthly simulation outputs so the two can never diverge
#   again.
#
# INPUTS  (one per SSP scenario)
#   outputs/<ssp>/simulation/ctmc_spark_monthly_2025_2075_rep1000.rds
#   where <ssp> in {ssp126, ssp245, ssp585}.
#
#   Required column: `p_month_major_mean` = monthly major-outbreak probability
#   p_ay,major (the Bernoulli component). Grouping keys: district_id, year, month.
#
# OUTPUTS
#   outputs/tables/tbl_ek_a4_poisson_bernoulli.csv
#   Console summary + a self-check on the Section 6.5.1 in-text numbers.
#
# DEFINITIONS  (per district x scenario, over all months of 2025-2075)
#   p_i             = p_month_major_mean_i  (clamped to [0,1], NA dropped)
#   P_exact         = 1 - prod(1 - p_i)                     # Eq. 12
#   Lambda          = sum(p_i)
#   P_poisson       = 1 - exp(-Lambda)                      # Eq. A.14
#   rel_diff_pct    = 100 * |P_exact - P_poisson| / P_exact
#   Sp2_over_2      = sum(p_i^2) / 2                         # Taylor 2nd-order bound
#   n_p_gt_001      = #{ p_i > 0.01 }
#   n_p_gt_005      = #{ p_i > 0.05 }
#   max_p           = max(p_i)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(here)
})

# ---- Configuration ----------------------------------------------------------
SSP_SCENARIOS <- c("ssp126", "ssp245", "ssp585")

SSP_LABELS <- c(
  ssp126 = "SSP1-2.6",
  ssp245 = "SSP2-4.5",
  ssp585 = "SSP5-8.5"
)

# Canonical district_id -> display label mapping (matches bulgular.Rmd / thesis).
DISTRICT_LABELS <- c(
  "TUR.40.25_1" = "Kartal (Istanbul)",
  "TUR.10.4_1"  = "Hopa (Artvin)",
  "TUR.59.4_1"  = "Zonguldak",
  "TUR.81.6_1"  = "Fethiye (Mugla)",
  "TUR.39.3_1"  = "Egirdir (Isparta)"
)

# Display order of districts within each scenario block (risk-descending, as in table).
DISTRICT_ORDER <- c("TUR.40.25_1", "TUR.10.4_1", "TUR.59.4_1", "TUR.81.6_1", "TUR.39.3_1")

MONTHLY_FILE <- "ctmc_spark_monthly_2025_2075_rep1000.rds"
PROB_COL     <- "p_month_major_mean"   # p_ay,major  (Bernoulli component)

# Resolve the monthly RDS path for a given scenario. Uses here() when available;
# falls back to a plain relative path if the project root is not detected.
monthly_path <- function(ssp) {
  rel <- file.path("outputs", ssp, "simulation", MONTHLY_FILE)
  p <- tryCatch(here::here(rel), error = function(e) rel)
  if (!file.exists(p)) p <- rel
  p
}

# ---- Core computation -------------------------------------------------------
# Given a data frame of monthly rows for ONE scenario, return one row per
# district with the Bernoulli-vs-Poisson comparison quantities.
compare_one_scenario <- function(monthly_df, ssp) {

  stopifnot(PROB_COL %in% names(monthly_df))
  stopifnot("district_id" %in% names(monthly_df))

  monthly_df %>%
    mutate(p = pmin(pmax(.data[[PROB_COL]], 0), 1)) %>%
    filter(!is.na(p)) %>%
    group_by(district_id) %>%
    summarise(
      n_months     = dplyr::n(),
      max_p        = max(p),
      n_p_gt_001   = sum(p > 0.01),
      n_p_gt_005   = sum(p > 0.05),
      P_exact      = 1 - prod(1 - p),          # Eq. 12  (heterogeneous Bernoulli)
      Lambda       = sum(p),
      Sp2_over_2   = sum(p^2) / 2,             # 2nd-order Taylor bound
      .groups = "drop"
    ) %>%
    mutate(
      P_poisson    = 1 - exp(-Lambda),         # Eq. A.14 (Poisson hazard)
      rel_diff_pct = dplyr::if_else(
        P_exact > 0,
        100 * abs(P_exact - P_poisson) / P_exact,
        NA_real_
      ),
      ssp          = ssp,
      ssp_label    = unname(SSP_LABELS[ssp]),
      district_label = dplyr::recode(district_id, !!!DISTRICT_LABELS)
    )
}

# ---- Load + compute across all scenarios ------------------------------------
message("Reading monthly simulation outputs and computing Table A.4 ...")

results <- purrr::map_dfr(SSP_SCENARIOS, function(ssp) {
  path <- monthly_path(ssp)
  if (!file.exists(path)) {
    warning(sprintf("[%s] file not found, skipping: %s", ssp, path))
    return(NULL)
  }
  df <- readRDS(path)
  message(sprintf("  [%s] %-55s  rows=%d", ssp, basename(path), nrow(df)))
  compare_one_scenario(df, ssp)
})

if (nrow(results) == 0) {
  stop("No scenario files could be read. Check the outputs/<ssp>/simulation/ paths.")
}

# ---- Order rows exactly as in the manuscript table --------------------------
results <- results %>%
  mutate(
    ssp          = factor(ssp, levels = SSP_SCENARIOS),
    district_id  = factor(district_id, levels = DISTRICT_ORDER)
  ) %>%
  arrange(ssp, district_id) %>%
  mutate(
    ssp         = as.character(ssp),
    district_id = as.character(district_id)
  )

# ---- Assemble the output table (column order matches Ek A.4) ----------------
tbl_a4 <- results %>%
  transmute(
    Senaryo        = ssp_label,
    Ilce           = district_label,
    `maks p_ay`    = max_p,
    `p>0,01 ay`    = n_p_gt_001,
    `p>0,05 ay`    = n_p_gt_005,
    `P_ufuk (kesin)`   = P_exact,
    `P_ufuk (Poisson)` = P_poisson,
    `Bagil fark (%)`   = rel_diff_pct,
    `Sp2/2`            = Sp2_over_2
  )

# ---- Write CSV --------------------------------------------------------------
out_dir <- tryCatch(here::here("outputs", "tables"), error = function(e) file.path("outputs", "tables"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "tbl_ek_a4_poisson_bernoulli.csv")
readr::write_csv(tbl_a4, out_csv)
message(sprintf("\nTable written to: %s", out_csv))

# ---- Console print (scientific, readable) -----------------------------------
fmt_sci <- function(x) formatC(x, format = "e", digits = 3)
fmt_pct <- function(x) formatC(x, format = "f", digits = 4)

cat("\n================ APPENDIX TABLE A.4 (reproduced) ================\n")
print_tbl <- results %>%
  transmute(
    Senaryo = ssp_label,
    Ilce    = district_label,
    max_p   = fmt_sci(max_p),
    `p>.01` = n_p_gt_001,
    `p>.05` = n_p_gt_005,
    P_exact = fmt_sci(P_exact),
    P_pois  = fmt_sci(P_poisson),
    `rel%`  = fmt_pct(rel_diff_pct),
    `Sp2/2` = fmt_sci(Sp2_over_2)
  )
print(as.data.frame(print_tbl), row.names = FALSE)

# ---- Summary of the maximum relative difference -----------------------------
worst <- results %>% arrange(desc(rel_diff_pct)) %>% slice(1)
cat(sprintf(
  "\nMaksimum bagil fark: %%%.3f  (%s, %s)\n",
  worst$rel_diff_pct, worst$district_label, worst$ssp_label
))

# ---- Self-check against Section 6.5.1 in-text numbers -----------------------
# The corrected prose should read (Kartal SSP5-8.5):
#   P_exact = 0.494, P_poisson = 0.4925, rel = 0.309%.
# Flag any drift so the prose is never again out of sync with the table.
cat("\n---- Bolum 6.5.1 metni ozdenetimi (Kartal SSP5-8.5) ----\n")
k585 <- results %>% filter(district_id == "TUR.40.25_1", ssp == "ssp585")
if (nrow(k585) == 1) {
  cat(sprintf("  P_ufuk (kesin)   = %.4f   -> metinde yazilmasi gereken: %.3f\n",
              k585$P_exact, round(k585$P_exact, 3)))
  cat(sprintf("  P_ufuk (Poisson) = %.4f   -> metinde yazilmasi gereken: %.4f\n",
              k585$P_poisson, round(k585$P_poisson, 4)))
  cat(sprintf("  Bagil fark       = %%%.3f\n", k585$rel_diff_pct))

  STALE <- c(exact = 0.484, pois = 0.482)   # values currently (wrongly) in prose
  if (abs(k585$P_exact - STALE["exact"]) > 1e-3 ||
      abs(k585$P_poisson - STALE["pois"]) > 1e-3) {
    cat("  [UYARI] Bolum 6.5.1 metnindeki 0,484 / 0,482 degerleri BAYAT.\n")
    cat("          Dogru degerlerle guncelleyin (yukariya bakiniz).\n")
  }
} else {
  cat("  [BILGI] Kartal SSP5-8.5 satiri bulunamadi (ssp585 dosyasi eksik olabilir).\n")
}

cat("\nBitti.\n")
