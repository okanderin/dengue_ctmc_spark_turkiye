# ==========================================================
# validation_two_stage.R
# Two-stage CTMC spark model validation
#
# Stage 1 — Software verification:
#   use_stochastic_EIP = FALSE (deterministic EIP mean)
#   MC P_est vs analytical gambler's ruin P_est
#   Expected difference: ~10⁻⁷ (MC sampling noise only)
#
# Stage 2 — Jensen bias assessment:
#   use_stochastic_EIP = TRUE vs FALSE
#   Measures E[P_est(EIP_i)] - P_est(E[EIP])
#   Expected: |ΔP_est| < 5% (per Appendix 1.3)
#
# SAFE: Does NOT modify any pipeline outputs.
#       Writes results to outputs/validation/ directory.
#
# Usage:
#   Sys.setenv(SSP_SCENARIO = "ssp245")
#   source("R/validation_two_stage.R")
#
# To run all three SSPs:
#   for (s in c("ssp126", "ssp245", "ssp585")) {
#     Sys.setenv(SSP_SCENARIO = s)
#     source("R/validation_two_stage.R")
#   }
# ==========================================================

cat("\n========================================\n")
cat("Two-stage CTMC Validation\n")
cat("========================================\n\n")

## ---- 0) Setup ----
source("R/01_setup/init.R")
source("R/03_models/parameter_functions.R")

# Source ctmc_spark.R to load all helper functions.
# Problem: its entrypoint runs run_ctmc_spark() when in globalenv().
# Solution: Pre-define a no-op version, source the file (entrypoint fires
# the no-op), then the real run_ctmc_spark overwrites it right after.
# But wait — source() defines run_ctmc_spark BEFORE the entrypoint block,
# so the entrypoint calls the REAL function. We need another approach.
#
# SAFEST approach: source with a wrapper environment, then copy functions out.
# The entrypoint guard `identical(environment(), globalenv())` will be FALSE
# in .val_env, so the pipeline won't run.
#
# BUT: ctmc_spark.R lines 53-64 check for fp_import and may stop() if the
# legacy file doesn't exist. We pre-set fp_import in .val_env to avoid this.

.val_env <- new.env(parent = globalenv())

# Pre-set fp_import so ctmc_spark.R's guard doesn't stop()
.val_env$fp_import <- file.path(DIR_PROCESSED_SSP,
                                "importation_pressure_monthly_2025_2075.rds")

# Suppress the source() calls inside ctmc_spark.R (already loaded above)
# by pre-setting a flag that init.R checks, or just let them re-run (harmless).
source("R/03_models/ctmc_spark.R", local = .val_env)

# Copy needed functions to global environment
for (fn_name in c("compute_ctmc_month_spark", "extinction_prob_bd",
                   "load_importation_monthly", "load_climate_monthly_long",
                   "load_or_build_species_map", "read_trait_params",
                   "read_rds_checked", "assert_has_cols", "assert_unique_key",
                   "assert_numeric_nonneg", "make_lambda_local",
                   ".major_rule_validate", ".sim_bd_once", ".sim_bd_prob",
                   "%||%")) {
  if (exists(fn_name, envir = .val_env, inherits = FALSE)) {
    assign(fn_name, get(fn_name, envir = .val_env), envir = globalenv())
  }
}
rm(.val_env)

cat("✓ Functions loaded: compute_ctmc_month_spark, extinction_prob_bd, etc.\n\n")

## ---- 1) Model parameters (must match main pipeline) ----
YEAR_MIN <- 2025L
YEAR_MAX <- 2075L
TAU      <- 30L
IP_DAYS  <- 5
BETA_VH  <- 0.3
BETA_HV  <- 0.33
M_RATIO  <- 1.0
EIP_MC_N <- 2000L
EIP_SEED <- 123L
GAMMA_DAY <- 1 / IP_DAYS

## ---- 2) Load data ----
fp_import_file  <- file.path(DIR_PROCESSED_SSP,
                             "importation_pressure_monthly_2025_2075.rds")
fp_climate_file <- file.path(DIR_PROCESSED_SSP,
                             "climate_sentinel_monthly_long_2015_2100.rds")
fp_traits_aeg   <- file.path(DIR_PROCESSED, "trait_params_aegypti.csv")
fp_traits_alb   <- file.path(DIR_PROCESSED, "trait_params_albopictus.csv")

# Check all inputs exist
for (fp in c(fp_import_file, fp_climate_file, fp_traits_aeg, fp_traits_alb)) {
  if (!file.exists(fp)) stop("Missing input: ", fp, call. = FALSE)
}

df_import  <- load_importation_monthly(fp_import_file)
df_climate <- load_climate_monthly_long(fp_climate_file)

df_import  <- df_import  %>% dplyr::filter(year >= YEAR_MIN, year <= YEAR_MAX)
df_climate <- df_climate %>% dplyr::filter(year >= YEAR_MIN, year <= YEAR_MAX)

key <- c("district_id", "year", "month")
df <- df_climate %>%
  dplyr::left_join(df_import, by = key)
df$lambda_import[is.na(df$lambda_import)]     <- 0
df$q_import_month[is.na(df$q_import_month)]   <- 0

sp_map <- load_or_build_species_map(unique(df$district_id))
df <- df %>% dplyr::left_join(sp_map, by = "district_id")

tab_aeg <- read_trait_params(fp_traits_aeg)
tab_alb <- read_trait_params(fp_traits_alb)

cat("Data loaded. Rows:", nrow(df), "\n")
cat("SSP scenario:", SSP_SCENARIO, "\n\n")

## ---- 3) Compute P_est for both configurations ----
cat("Computing P_est (deterministic EIP & stochastic EIP)...\n")
cat("This may take 5-15 minutes per SSP.\n\n")

n_rows <- nrow(df)
p_est_det   <- numeric(n_rows)   # use_stochastic_EIP = FALSE
p_est_stoch <- numeric(n_rows)   # use_stochastic_EIP = TRUE
lam_det     <- numeric(n_rows)   # lambda_local_i1 (deterministic)
lam_stoch   <- numeric(n_rows)   # mean lambda_local_i1 (stochastic)

t0 <- Sys.time()

for (i in seq_len(n_rows)) {
  row <- df[i, ]
  tab_traits <- if (row$species == "aegypti") tab_aeg else tab_alb
  
  # Deterministic seed (same as main pipeline)
  key_str   <- paste(row$district_id, row$year, row$month, sep = "|")
  seed_here <- as.integer((abs(sum(utf8ToInt(key_str))) + EIP_SEED) %% .Machine$integer.max)
  if (seed_here == 0L) seed_here <- 1L
  
  # --- Stage 1: Deterministic EIP (σ_EIP = 0) ---
  res_det <- compute_ctmc_month_spark(
    temp_c = row$temp_c, rh = row$rh,
    lambda_import_month = row$lambda_import,
    tab_traits = tab_traits,
    beta_vh = BETA_VH, beta_hv = BETA_HV, m = M_RATIO,
    gamma_day = GAMMA_DAY, tau = TAU, i0 = 1L,
    include_import_in_q = FALSE,
    use_stochastic_EIP = FALSE,
    eip_mc_n = EIP_MC_N, eip_seed = seed_here
  )
  p_est_det[i] <- res_det$p_establishment
  lam_det[i]   <- res_det$lambda_local_i1
  
  # --- Stage 2: Stochastic EIP (σ_log ≈ 0.45) ---
  res_stoch <- compute_ctmc_month_spark(
    temp_c = row$temp_c, rh = row$rh,
    lambda_import_month = row$lambda_import,
    tab_traits = tab_traits,
    beta_vh = BETA_VH, beta_hv = BETA_HV, m = M_RATIO,
    gamma_day = GAMMA_DAY, tau = TAU, i0 = 1L,
    include_import_in_q = FALSE,
    use_stochastic_EIP = TRUE,
    eip_mc_n = EIP_MC_N, eip_seed = seed_here
  )
  p_est_stoch[i] <- res_stoch$p_establishment
  lam_stoch[i]   <- res_stoch$lambda_local_i1
  
  # Progress
  if (i %% 500 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    pct <- round(i / n_rows * 100, 1)
    eta <- round(elapsed / i * (n_rows - i), 1)
    cat(sprintf("  [%d/%d] %s%%  elapsed=%.1f min  ETA=%.1f min\n",
                i, n_rows, pct, elapsed, eta))
  }
}

elapsed_total <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nDone. Total time: %.1f minutes\n\n", elapsed_total))

## ---- 4) Analytical reference (gambler's ruin formula) ----
# P_est_analytical = (1 - rho) / (1 - rho^tau)  where rho = gamma / lambda
# This uses the DETERMINISTIC lambda (no EIP stochasticity)
p_est_analytical <- numeric(n_rows)

for (i in seq_len(n_rows)) {
  li <- lam_det[i]
  if (li <= 0) {
    p_est_analytical[i] <- 0
    next
  }
  rho <- GAMMA_DAY / li
  if (abs(rho - 1) < 1e-12) {
    p_est_analytical[i] <- 1 / TAU
  } else {
    p_est_analytical[i] <- (1 - rho) / (1 - rho^TAU)
  }
  # Clamp to [0, 1]
  p_est_analytical[i] <- max(0, min(1, p_est_analytical[i]))
}

## ---- 5) Assemble results ----
df_val <- df %>%
  dplyr::mutate(
    p_est_det        = p_est_det,
    p_est_stoch      = p_est_stoch,
    p_est_analytical = p_est_analytical,
    lambda_det       = lam_det,
    lambda_stoch     = lam_stoch,
    
    # Stage 1: MC (det EIP) vs Analytical
    diff_verify      = p_est_det - p_est_analytical,
    abs_diff_verify  = abs(p_est_det - p_est_analytical),
    
    # Stage 2: Stochastic EIP vs Deterministic EIP (Jensen bias)
    diff_jensen      = p_est_stoch - p_est_det,
    abs_diff_jensen  = abs(p_est_stoch - p_est_det),
    pct_diff_jensen  = dplyr::if_else(
      p_est_det > 1e-15,
      (p_est_stoch - p_est_det) / p_est_det * 100,
      NA_real_
    )
  )

## ---- 6) Stage 1 Summary: Software Verification ----
active_rows <- df_val %>% dplyr::filter(lambda_det > 0)

tbl_stage1 <- active_rows %>%
  dplyr::summarise(
    ssp           = SSP_SCENARIO,
    n_total       = nrow(df_val),
    n_active      = dplyr::n(),
    mean_abs_diff = mean(abs_diff_verify, na.rm = TRUE),
    max_abs_diff  = max(abs_diff_verify, na.rm = TRUE),
    pass_002      = sum(abs_diff_verify < 0.02, na.rm = TRUE),
    pass_pct      = round(sum(abs_diff_verify < 0.02, na.rm = TRUE) / dplyr::n() * 100, 1)
  )

cat("=== STAGE 1: Software Verification (σ_EIP = 0) ===\n")
cat("MC P_est(deterministic EIP) vs Analytical Gambler's Ruin\n\n")
print(as.data.frame(tbl_stage1))
cat("\n")

## ---- 7) Stage 2 Summary: Jensen Bias Assessment ----
tbl_stage2 <- active_rows %>%
  dplyr::filter(p_est_det > 1e-15) %>%   # exclude near-zero for % calc

  dplyr::summarise(
    ssp              = SSP_SCENARIO,
    n_active         = dplyr::n(),
    mean_abs_diff    = mean(abs_diff_jensen, na.rm = TRUE),
    max_abs_diff     = max(abs_diff_jensen, na.rm = TRUE),
    mean_pct_jensen  = mean(pct_diff_jensen, na.rm = TRUE),
    median_pct_jensen = stats::median(pct_diff_jensen, na.rm = TRUE),
    max_abs_pct_jensen = max(abs(pct_diff_jensen), na.rm = TRUE),
    within_5pct      = sum(abs(pct_diff_jensen) < 5, na.rm = TRUE),
    within_5pct_pct  = round(sum(abs(pct_diff_jensen) < 5, na.rm = TRUE) / dplyr::n() * 100, 1)
  )

cat("=== STAGE 2: Jensen Bias Assessment (σ_log ≈ 0.45) ===\n")
cat("E[P_est(EIP_i)] vs P_est(E[EIP])\n\n")
print(as.data.frame(tbl_stage2))
cat("\n")

## ---- 8) Save outputs ----
out_dir <- file.path(DIR_OUTPUT_SSP, "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Full row-level results
fp_full <- file.path(out_dir, "validation_two_stage_full.rds")
saveRDS(df_val, fp_full)
cat("Full results saved:", fp_full, "\n")

# Summary tables
fp_s1 <- file.path(out_dir, "validation_stage1_summary.csv")
fp_s2 <- file.path(out_dir, "validation_stage2_jensen.csv")
readr::write_csv(tbl_stage1, fp_s1)
readr::write_csv(tbl_stage2, fp_s2)
cat("Stage 1 summary:", fp_s1, "\n")
cat("Stage 2 summary:", fp_s2, "\n")

# Row-level Jensen bias (for plotting)
fp_jensen_detail <- file.path(out_dir, "validation_jensen_detail.csv")
jensen_detail <- df_val %>%
  dplyr::filter(lambda_det > 0, p_est_det > 1e-15) %>%
  dplyr::select(district_id, year, month, temp_c, rh,
                p_est_det, p_est_stoch, p_est_analytical,
                diff_verify, diff_jensen, pct_diff_jensen)
readr::write_csv(jensen_detail, fp_jensen_detail)
cat("Jensen detail:", fp_jensen_detail, "\n")

cat("\n========================================\n")
cat("Validation complete for", SSP_SCENARIO, "\n")
cat("========================================\n")
