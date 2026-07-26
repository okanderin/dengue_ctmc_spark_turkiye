# =========================================================
# R/04_results/6.10_pest_regression.R   (CORRECTED)
# ---------------------------------------------------------
# T_opt estimation for P_est, fixing three prior issues:
#
#   (K5a) Selection on the dependent variable.
#         The old script filtered `P_est > 1e-15` where the DV is
#         log10(P_est) -> truncation bias that mechanically shifts T_opt.
#         FIX: keep ALL district-months; treat P_est <= floor as
#         LEFT-CENSORED at log10(floor) and fit a Tobit (survreg).
#
#   (K5c) Random slope on 5 clusters.
#         (T | district_id) with only 5 districts gives unreliable
#         variance components (>=30-cluster heuristic; boundary chi^2).
#         FIX: district FIXED effects are primary. Random-slope model
#         removed. A random-INTERCEPT LMM is kept ONLY as a caveated
#         secondary block for ICC/R^2 continuity.
#
#   (K5b) Circular "validation" against Mordecai.
#         Mordecai's thermal curves are the model INPUT, so agreement
#         verifies coding, not external validity.
#         FIX: PRIMARY T_opt is now MECHANISTIC (argmax_T of lambda(i=1)),
#         and the Mordecai comparison is relabelled "internal consistency
#         check". A real hindcast validation belongs elsewhere.
#
# Output (under outputs/tables/):
#   topt_mechanistic.csv        (PRIMARY, bias-free, per district + species)
#   tobit_fe_coefficients.csv   (censored, all months, district FE)
#   ols_active_diagnostic.csv   (active-only OLS-FE = truncation-bias probe)
#   lmm_ri_secondary.csv        (random-intercept ICC/R^2, 5-cluster caveat)
#   cross_ssp_topt.csv          (mechanistic + Tobit across 3 SSPs)
#   topt_comparison.csv         (internal-consistency summary)
# =========================================================

source("R/04_results/00_results_setup.R")

suppressPackageStartupMessages({
  library(survival)     # survreg() Tobit (base recommended pkg)
  library(lme4)
  library(lmerTest)
  library(performance)  # ICC + Nakagawa R^2
  library(broom)
  library(broom.mixed)
})

# ---- Settings ----
TARGET_SSP  <- "ssp245"
ALL_SSPS    <- c("ssp126", "ssp245", "ssp585")
P_EST_FLOOR <- 1e-15
LOG_FLOOR   <- log10(P_EST_FLOOR)          # = -15
TABLES_DIR  <- here("outputs", "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

# Baseline biology (constants — do NOT affect argmax_T of lambda, only its scale)
BETA_VH <- 0.30
BETA_HV <- 0.33
M_BASE  <- 1.00

# =========================================================
# Helper: load FULL monthly panel (NO truncation filter)
#   Non-finite / zero P_est kept and marked left-censored at the floor.
# =========================================================
load_pest_panel <- function(ssp) {
  d       <- load_ssp_data(ssp)
  monthly <- d$monthly

  p_col <- if ("p_establishment_mean" %in% names(monthly)) {
    "p_establishment_mean"
  } else if ("p_establishment" %in% names(monthly)) {
    "p_establishment"
  } else {
    stop("No P_est column in monthly output for ", ssp,
         "\nAvailable: ", paste(names(monthly), collapse = ", "))
  }

  required <- c("district_id", "year", "month", p_col, "temp_c", "rh")
  missing  <- setdiff(required, names(monthly))
  if (length(missing) > 0) {
    stop("Missing columns for ", ssp, ": ", paste(missing, collapse = ", "))
  }

  monthly %>%
    rename(P_est = !!p_col, T = temp_c, RH = rh) %>%
    mutate(
      ssp        = ssp,
      P_est      = ifelse(is.finite(P_est), P_est, 0),
      censored   = P_est <= P_EST_FLOOR,                 # left-censored flag
      log10_Pest = log10(pmax(P_est, P_EST_FLOOR)),      # censored rows sit at LOG_FLOOR
      event      = as.integer(!censored)                 # survreg: 1=observed, 0=left-censored
    )
}

topt_from_quad <- function(coef_vec) {
  bT  <- coef_vec[["T"]]
  bT2 <- coef_vec[["I(T^2)"]]
  as.numeric(-bT / (2 * bT2))
}

# =========================================================
# (1) MECHANISTIC T_opt  —  PRIMARY, bias-free
#     P_est is monotone increasing in lambda(i=1); therefore
#     argmax_T P_est = argmax_T lambda(i=1, T, RH).
#     m / beta are multiplicative constants -> do NOT move the peak.
# =========================================================
cat("\n========== (1) MECHANISTIC T_opt (argmax lambda) ==========\n")

mech_topt <- tryCatch({
  source("R/03_models/parameter_functions.R")  # read_trait_params(), make_lambda_local()

  tab_aeg <- read_trait_params(here("data_processed", "trait_params_aegypti.csv"))
  tab_alb <- read_trait_params(here("data_processed", "trait_params_albopictus.csv"))

  sp_path <- here("data_processed", "sentinel_species.csv")
  species_map <- if (file.exists(sp_path)) {
    readr::read_csv(sp_path, show_col_types = FALSE) %>%
      select(district_id, species)
  } else {
    tibble(district_id = character(), species = character())
  }

  lambda_curve <- function(T_vec, RH, tab) {
    vapply(T_vec, function(t) {
      make_lambda_local(n_infected = 1L, T = t, RH = RH,
                        beta_vh = BETA_VH, beta_hv = BETA_HV, m = M_BASE,
                        tab_traits = tab, use_stochastic_EIP = FALSE)
    }, numeric(1))
  }

  T_grid <- seq(10, 42, by = 0.02)

  # Representative RH per district = mean RH over the active season (TARGET_SSP)
  rh_by_dist <- load_pest_panel(TARGET_SSP) %>%
    filter(!censored) %>%
    group_by(district_id) %>%
    summarise(RH_rep = mean(RH, na.rm = TRUE), .groups = "drop")

  rh_by_dist %>%
    left_join(species_map, by = "district_id") %>%
    mutate(species = ifelse(is.na(species), "albopictus", species)) %>%
    rowwise() %>%
    mutate(
      T_opt_mech = {
        tab <- if (species == "aegypti") tab_aeg else tab_alb
        lam <- lambda_curve(T_grid, RH_rep, tab)
        if (all(!is.finite(lam)) || max(lam, na.rm = TRUE) <= 0) NA_real_
        else T_grid[which.max(lam)]
      }
    ) %>%
    ungroup()
}, error = function(e) {
  message("⚠ Mekanistik T_opt hesaplanamadi (parameter_functions / trait dosyalari): ",
          conditionMessage(e))
  NULL
})

T_opt_mech_overall <- NA_real_
if (!is.null(mech_topt)) {
  print(as.data.frame(mech_topt))
  T_opt_mech_overall <- round(mean(mech_topt$T_opt_mech, na.rm = TRUE), 2)
  cat(sprintf("Mekanistik T_opt (ilce ortalamasi) = %.2f °C\n", T_opt_mech_overall))
  write_csv(mech_topt, file.path(TABLES_DIR, "topt_mechanistic.csv"))
} else {
  cat("Mekanistik T_opt atlandi; regresyon T_opt'lari yine de raporlanacak.\n")
}

# =========================================================
# (2) DESCRIPTIVE regression — TOBIT, all months, district FIXED effects
#     Left-censored at LOG_FLOOR. Removes the truncation bias (K5a) and
#     uses fixed effects instead of 5-cluster random effects (K5c).
# =========================================================
cat("\n========== (2) Tobit-FE (censored, all months) ==========\n")

panel <- load_pest_panel(TARGET_SSP)
cat("Panel n =", nrow(panel),
    " | sansurlu (P_est<=floor) =", sum(panel$censored),
    " | gozlenen =", sum(!panel$censored),
    " | ilce =", n_distinct(panel$district_id), "\n")

tobit_fit <- survreg(
  Surv(log10_Pest, event, type = "left") ~ T + I(T^2) + RH + factor(district_id),
  data = panel, dist = "gaussian"
)

T_opt_tobit <- topt_from_quad(coef(tobit_fit))
cat(sprintf("Tobit-FE: beta_T = %.3f, beta_T2 = %.4f, T_opt = %.2f °C\n",
            coef(tobit_fit)[["T"]], coef(tobit_fit)[["I(T^2)"]], T_opt_tobit))

tobit_tab <- broom::tidy(tobit_fit) %>%
  filter(term %in% c("(Intercept)", "T", "I(T^2)", "RH")) %>%
  transmute(term,
            estimate = round(estimate, 4),
            std.error = round(std.error, 4),
            statistic = round(statistic, 3),
            p.value  = format.pval(p.value, digits = 3, eps = 1e-4))
print(tobit_tab)
write_csv(tobit_tab, file.path(TABLES_DIR, "tobit_fe_coefficients.csv"))

# =========================================================
# (3) DIAGNOSTIC — active-only OLS-FE (quantifies truncation bias)
#     Same DV/covariates but ONLY observed months, to expose the gap
#     that the old (filtered) analysis suffered.
# =========================================================
cat("\n========== (3) OLS-FE active-only (truncation-bias probe) ==========\n")

ols_active <- lm(log10_Pest ~ T + I(T^2) + RH + factor(district_id),
                 data = subset(panel, !censored))
T_opt_ols_active <- topt_from_quad(coef(ols_active))
bias_deg <- T_opt_ols_active - T_opt_tobit

cat(sprintf("OLS-FE (yalnizca aktif ay): T_opt = %.2f °C  (n=%d)\n",
            T_opt_ols_active, sum(!panel$censored)))
cat(sprintf("Truncation bias = T_opt(aktif) - T_opt(Tobit) = %.2f °C\n", bias_deg))

ols_active_tab <- broom::tidy(ols_active, conf.int = TRUE) %>%
  filter(term %in% c("(Intercept)", "T", "I(T^2)", "RH")) %>%
  transmute(term,
            estimate = round(estimate, 4),
            std.error = round(std.error, 4),
            conf.low = round(conf.low, 4),
            conf.high = round(conf.high, 4),
            p.value  = format.pval(p.value, digits = 3, eps = 1e-4)) %>%
  mutate(R2_active = round(summary(ols_active)$r.squared, 3),
         T_opt_active = round(T_opt_ols_active, 2),
         truncation_bias_deg = round(bias_deg, 2))
print(ols_active_tab)
write_csv(ols_active_tab, file.path(TABLES_DIR, "ols_active_diagnostic.csv"))

# =========================================================
# (4) SECONDARY — random-INTERCEPT LMM (ICC / Nakagawa R^2)
#     Kept ONLY for continuity with prior tables. NO random slope.
#     CAVEAT: variance components rest on 5 clusters (< 30) and are
#     approximate; NOT used for the primary T_opt.
# =========================================================
cat("\n========== (4) Random-intercept LMM (SECONDARY, 5-cluster caveat) ==========\n")

lmm_ri <- lmer(log10_Pest ~ T + I(T^2) + RH + (1 | district_id),
               data = subset(panel, !censored), REML = TRUE)

icc_val <- performance::icc(lmm_ri)
r2_val  <- performance::r2_nakagawa(lmm_ri)
cat("ICC (adjusted):", round(icc_val$ICC_adjusted, 3),
    " | R2_marg:", round(r2_val$R2_marginal, 3),
    " | R2_cond:", round(r2_val$R2_conditional, 3), "\n")
cat("UYARI: 5 ilce (kume) uzerinden; varyans bilesenleri yaklasiktir (>=30 kume onerilir).\n")

lmm_ri_tab <- tibble(
  metric = c("ICC_adjusted", "R2_marginal", "R2_conditional",
             "n_clusters", "note"),
  value  = c(round(icc_val$ICC_adjusted, 3),
             round(r2_val$R2_marginal, 3),
             round(r2_val$R2_conditional, 3),
             n_distinct(subset(panel, !censored)$district_id),
             NA_real_)
) %>%
  mutate(note = c(rep("", 4),
                  "5 kume: varyans bilesenleri yaklasik; birincil T_opt DEGIL")[seq_len(n())])
print(lmm_ri_tab)
write_csv(lmm_ri_tab, file.path(TABLES_DIR, "lmm_ri_secondary.csv"))

# =========================================================
# (5) Cross-SSP T_opt  —  mechanistic + Tobit-FE across scenarios
# =========================================================
cat("\n========== (5) Cross-SSP T_opt ==========\n")

cross_ssp <- lapply(ALL_SSPS, function(s) {
  p  <- load_pest_panel(s)
  tf <- survreg(Surv(log10_Pest, event, type = "left") ~ T + I(T^2) + RH +
                  factor(district_id), data = p, dist = "gaussian")
  # mechanistic (dominant species + overall mean RH) if available
  t_mech <- NA_real_
  if (!is.null(mech_topt) && exists("make_lambda_local")) {
    rh_bar <- mean(p$RH[!p$censored], na.rm = TRUE)
    tabo   <- tryCatch(read_trait_params(here("data_processed",
                                              "trait_params_albopictus.csv")),
                       error = function(e) NULL)
    if (!is.null(tabo)) {
      Tg  <- seq(10, 42, by = 0.02)
      lam <- vapply(Tg, function(t)
        make_lambda_local(1L, t, rh_bar, beta_vh = BETA_VH, beta_hv = BETA_HV,
                          m = M_BASE, tab_traits = tabo,
                          use_stochastic_EIP = FALSE), numeric(1))
      if (max(lam, na.rm = TRUE) > 0) t_mech <- Tg[which.max(lam)]
    }
  }
  tibble(
    SSP_code       = s,
    SSP_label      = SSP_LABELS[s],
    n_total        = nrow(p),
    n_censored     = sum(p$censored),
    T_opt_mech_alb = round(t_mech, 2),
    T_opt_tobit    = round(topt_from_quad(coef(tf)), 2)
  )
}) %>% bind_rows()
print(as.data.frame(cross_ssp))
write_csv(cross_ssp, file.path(TABLES_DIR, "cross_ssp_topt.csv"))

# =========================================================
# (6) Comparison / internal-consistency summary
#     Mordecai comparison relabelled: NOT validation.
# =========================================================
cat("\n========== (6) T_opt karsilastirmasi ==========\n")

topt_comparison <- tibble(
  Yontem = c(
    "Mekanistik (∂λ/∂T=0) — BIRINCIL, yansiz",
    "Tobit-FE (sansurlu, tum aylar, ilce sabit etkisi)",
    "OLS-FE (yalnizca aktif ay) — TRUNCATION etkili (tani)",
    "Truncation bias (aktif − Tobit)",
    "Mordecai ~26–28°C ile karsilastirma"
  ),
  Deger = c(
    ifelse(is.na(T_opt_mech_overall), "—", sprintf("%.2f °C", T_opt_mech_overall)),
    sprintf("%.2f °C", T_opt_tobit),
    sprintf("%.2f °C", T_opt_ols_active),
    sprintf("%+.2f °C", bias_deg),
    "IC TUTARLILIK KONTROLU"
  ),
  Not = c(
    "P_est ithalattan bagimsiz; peak = argmax lambda; m/beta peaki kaydirmaz",
    "Sansur sinir alti aylar bilgi tasir; yansiz T_opt",
    "Eski filtrenin (P_est>1e-15) urettigi yanliligi gosterir",
    "Bu fark ne kadar buyukse eski T_opt o kadar yanliydi",
    "Mordecai egrileri modelin GIRDISI; uyum KODLAMA dogrulugunu gosterir, DIS gecerliligi DEGIL. Gercek dogrulama: Avrupa otokton salgin hindcast'i ([164] Farooq)."
  )
)
print(as.data.frame(topt_comparison))
write_csv(topt_comparison, file.path(TABLES_DIR, "topt_comparison.csv"))

cat("\n========== Bitti. Tablolar:", TABLES_DIR, "==========\n")
