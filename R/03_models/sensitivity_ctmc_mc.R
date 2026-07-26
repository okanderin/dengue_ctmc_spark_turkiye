# ==========================================================
# R/03_models/sensitivity_ctmc_mc.R
# Sensitivity analysis for CTMC Spark model — Monte Carlo version
# SSP-AWARE: outputs to DIR_OUTPUT_SSP/sensitivity/ctmc_mc
#
# Base scenario: m = 1.0 (Ross-Macdonald standard)
# Varied parameters (one-at-a-time):
#   m          : mosquito/human ratio         (0.50, 0.80, 1.00, 1.20, 2.00)
#                0.50 = low-density urban, 0.80 = early colonisation,
#                1.00 = base (R-M standard), 1.20 = moderate, 2.00 = high-density
#   beta_vh    : vector→human transmission    (0.24, 0.30, 0.36)  [±20%]
#   ip_days    : infectious period [days]     (4, 5, 6)           [±20%]
#                (ip ↑ → γ ↓ → higher risk)
#
# MC settings (matched to main run):
#   use_stochastic_EIP = TRUE
#   n_mc_eip = 2000
#   n_rep    = 1000
#   seed     = 123
#
# Outputs (each scenario → own sub-folder):
#   outputs/sensitivity/ctmc_mc/<scenario>/ctmc_spark_*_rep1000.csv/.rds
#   outputs/sensitivity/ctmc_mc/sensitivity_summary.csv
#   outputs/sensitivity/ctmc_mc/sensitivity_tornado.csv
#
# Run order:
#   source("R/03_models/ctmc_spark_monte_carlo.R")  # defines run_ctmc_spark()
#   source("R/03_models/sensitivity_ctmc_mc.R")      # this file
# ==========================================================

run_sensitivity_mc <- function(
    year_min    = 2025,
    year_max    = 2075,
    n_rep       = 1000L,
    n_mc_eip    = 2000L,
    seed        = 123L,
    out_base    = file.path(DIR_OUTPUT_SSP, "sensitivity", "ctmc_mc")
) {
  
  if (!dir.exists(out_base)) dir.create(out_base, recursive = TRUE)
  
  # ----------------------------------------------------------
  # Base values — MUST match ctmc_spark_monte_carlo.R entrypoint
  # ----------------------------------------------------------
  BASE <- list(
    m       = 1.0,    # Conditional limit
    beta_vh = 0.3,    # vector-to-human transmission probability
    beta_hv = 0.33,   # human-to-vector transmission probability (fixed)
    ip_days = 5       # infectious period days; gamma = 1/ip_days
  )
  
  # Fixed parameters (not varied)
  TAU              <- 30L
  USE_STOCH_EIP    <- TRUE
  
  # ----------------------------------------------------------
  # Scenario grid: one-at-a-time (m: 0.50/0.80/1.00/1.20/2.00; others: ±20%)
  # ----------------------------------------------------------
  # ip_days note:
  #   ip_plus_20  = 6 days  →  gamma = 1/6 = 0.167/day  →  LOWER recovery → HIGHER risk
  #   ip_minus_20 = 4 days  →  gamma = 1/4 = 0.250/day  →  HIGHER recovery → LOWER risk
  
  # m: five biologically motivated values
  #   0.50 = low-density / sparse urban vector population
  #   0.80 = early colonisation / suburban fringe
  #   1.00 = base (Ross-Macdonald standard, 1 female/person)
  #   1.20 = moderate urban density (established Aedes population)
  #   2.00 = high-density tropical urban (upper plausible range)
  # beta_vh and ip_days: ±20% one-at-a-time
  scenarios <- list(
    base          = list(m = 1.00,   beta_vh = BASE$beta_vh,       ip_days = BASE$ip_days),
    m_050         = list(m = 0.50,   beta_vh = BASE$beta_vh,       ip_days = BASE$ip_days),
    m_080         = list(m = 0.80,   beta_vh = BASE$beta_vh,       ip_days = BASE$ip_days),
    m_120         = list(m = 1.20,   beta_vh = BASE$beta_vh,       ip_days = BASE$ip_days),
    m_200         = list(m = 2.00,   beta_vh = BASE$beta_vh,       ip_days = BASE$ip_days),
    beta_minus_20 = list(m = BASE$m, beta_vh = BASE$beta_vh * 0.8, ip_days = BASE$ip_days),
    beta_plus_20  = list(m = BASE$m, beta_vh = BASE$beta_vh * 1.2, ip_days = BASE$ip_days),
    ip_plus_20    = list(m = BASE$m, beta_vh = BASE$beta_vh,       ip_days = BASE$ip_days * 1.2),
    ip_minus_20   = list(m = BASE$m, beta_vh = BASE$beta_vh,       ip_days = BASE$ip_days * 0.8)
  )
  
  # ----------------------------------------------------------
  # Parameter table (for console + CSV)
  # ----------------------------------------------------------
  param_tbl <- purrr::imap_dfr(scenarios, function(s, name) {
    tibble::tibble(
      scenario  = name,
      m         = s$m,
      beta_vh   = s$beta_vh,
      beta_hv   = BASE$beta_hv,
      ip_days   = s$ip_days,
      gamma_day = round(1 / s$ip_days, 4)
    )
  })
  
  message("\n=== Sensitivity scenario parameters ===")
  print(param_tbl, n = Inf)
  readr::write_csv(param_tbl, file.path(out_base, "sensitivity_params.csv"))
  
  # ----------------------------------------------------------
  # Run each scenario
  # ----------------------------------------------------------
  all_horizon <- list()
  
  for (name in names(scenarios)) {
    s <- scenarios[[name]]
    message("\n", strrep("=", 55))
    message("Scenario: ", name,
            "  (m=", s$m,
            ", beta_vh=", s$beta_vh,
            ", ip_days=", s$ip_days, ")")
    message(strrep("=", 55))
    
    # Each scenario writes to its own sub-directory
    # → prevents file overwriting between scenarios
    scen_dir <- file.path(out_base, name)
    if (!dir.exists(scen_dir)) dir.create(scen_dir, recursive = TRUE)
    
    res <- run_ctmc_spark(
      year_min               = year_min,
      year_max               = year_max,
      tau                    = TAU,
      infectious_period_days = s$ip_days,
      beta_vh                = s$beta_vh,
      beta_hv                = BASE$beta_hv,
      m                      = s$m,
      include_import_in_q    = FALSE,
      use_stochastic_EIP     = USE_STOCH_EIP,
      n_mc_eip               = n_mc_eip,
      n_rep                  = n_rep,
      seed                   = seed,
      major_rule             = list(type = "establishment"),
      out_dir                = scen_dir
    )
    
    all_horizon[[name]] <- res$horizon %>%
      dplyr::mutate(
        scenario  = name,
        m_val     = s$m,
        beta_vh   = s$beta_vh,
        beta_hv   = BASE$beta_hv,
        ip_days   = s$ip_days,
        gamma_day = 1 / s$ip_days
      )
  }
  
  # ----------------------------------------------------------
  # Combined summary
  # ----------------------------------------------------------
  summary_df <- dplyr::bind_rows(all_horizon)
  
  # Attach base values for % change calculation
  base_vals <- all_horizon[["base"]] %>%
    dplyr::select(district_id,
                  p_base_mean = p_ge1_major_mean,
                  p_base_p2_5 = p_ge1_major_p2_5,
                  p_base_p97_5 = p_ge1_major_p97_5)
  
  summary_df <- summary_df %>%
    dplyr::left_join(base_vals, by = "district_id") %>%
    dplyr::mutate(
      delta_abs = p_ge1_major_mean - p_base_mean,
      delta_pct = dplyr::if_else(
        p_base_mean > 1e-12,
        (p_ge1_major_mean - p_base_mean) / p_base_mean * 100,
        NA_real_
      )
    )
  
  readr::write_csv(summary_df, file.path(out_base, "sensitivity_summary.csv"))
  
  # ----------------------------------------------------------
  # Tornado input: mean % change across districts per scenario
  # ----------------------------------------------------------
  tornado <- summary_df %>%
    dplyr::filter(scenario != "base") %>%
    dplyr::group_by(scenario, m_val, beta_vh, ip_days, gamma_day) %>%
    dplyr::summarise(
      mean_delta_pct  = mean(delta_pct,        na.rm = TRUE),
      mean_delta_abs  = mean(delta_abs,         na.rm = TRUE),
      mean_p_mean     = mean(p_ge1_major_mean,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(abs(mean_delta_pct)))
  
  readr::write_csv(tornado, file.path(out_base, "sensitivity_tornado.csv"))
  
  message("\n", strrep("=", 55))
  message("SENSITIVITY ANALYSIS COMPLETE")
  message(strrep("=", 55))
  message("Outputs: ", out_base)
  
  message("\nTornado summary (mean % change vs base across districts):")
  print(tornado %>% dplyr::select(scenario, m_val, beta_vh, ip_days, mean_delta_pct), n = Inf)
  
  invisible(list(summary = summary_df, tornado = tornado))
}

# ==========================================================
# Entrypoint
# ==========================================================

# Step 1: Load the MC runner (defines run_ctmc_spark)
source("R/03_models/ctmc_spark_monte_carlo.R")

# Step 2: Run sensitivity analysis
# WARNING: n_rep=1000 x 9 scenarios is slow (~hours).
# For a quick test, set n_rep=50 first to verify it runs.
sensitivity_mc_res <- run_sensitivity_mc(
  year_min  = 2025,
  year_max  = 2075,
  n_rep     = 1000L,   # match main run
  n_mc_eip  = 2000L,   # match main run
  seed      = 123L,    # match main run
  out_base  = file.path(DIR_OUTPUT_SSP, "sensitivity", "ctmc_mc")
)