# =========================================================
# R/04_results/6.10_pest_regression.R
# ---------------------------------------------------------
# OLS (yearly) + LMM (monthly) regressions on log10(P_est)
#
# Methodological note:
#   Earlier analysis used log10(P_est) for OLS but log10(p_ay) for LMM, with
#   lambda_import as a regressor. Two issues:
#     (1) Different DV (P_est vs p_ay) prevents direct comparison.
#     (2) p_ay = 1 - exp(-lambda_import * P_est), so log10(p_ay) ~ log10(lambda_import)
#         is structurally tautological.
#
#   Both regressions now use log10(P_est). By gambler's-ruin construction,
#   P_est = (1 - rho)/(1 - rho^tau), rho = 1/R0, so P_est is a function of
#   R0(T, RH, biology) only and does not depend on lambda_import by definition.
#
# Inclusion rule:
#   Active observations only (P_est > 1e-15) -- avoids log10(0) = -Inf.
#
# Output (under outputs/tables/):
#   ols_hierarchical.csv      (Table 24)
#   lmm_fit_metrics.csv       (Table 25)
#   lmm_fixed_effects.csv     (Table 26)
#   lmm_r2_nakagawa.csv       (marginal/conditional R^2)
#   ols_lmm_comparison.csv    (Table 27)
#   cross_ssp_topt.csv        (Table 28)
# =========================================================

source("R/04_results/00_results_setup.R")

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)    # p-values for fixed effects
  library(performance) # ICC and Nakagawa R^2
  library(broom.mixed)
})

# ---- Settings ----
TARGET_SSP   <- "ssp245"
ALL_SSPS     <- c("ssp126", "ssp245", "ssp585")
P_EST_FLOOR  <- 1e-15
TABLES_DIR   <- here("outputs", "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers ----
load_pest_panel <- function(ssp) {
  d <- load_ssp_data(ssp)
  monthly <- d$monthly
  
  p_col <- if ("p_establishment_mean" %in% names(monthly)) {
    "p_establishment_mean"
  } else if ("p_establishment" %in% names(monthly)) {
    "p_establishment"
  } else {
    stop("No P_est column found in monthly output for ", ssp,
         "\nAvailable columns: ", paste(names(monthly), collapse = ", "))
  }
  
  required <- c("district_id", "year", "month", p_col, "temp_c", "rh")
  missing  <- setdiff(required, names(monthly))
  if (length(missing) > 0) {
    stop("Missing columns in monthly output for ", ssp, ": ",
         paste(missing, collapse = ", "),
         "\nAvailable: ", paste(names(monthly), collapse = ", "))
  }
  
  monthly %>%
    rename(P_est = !!p_col, T = temp_c, RH = rh) %>%
    mutate(ssp = ssp) %>%
    filter(P_est > P_EST_FLOOR, is.finite(P_est)) %>%
    mutate(log10_Pest = log10(P_est))
}

load_pest_yearly <- function(ssp) {
  load_pest_panel(ssp) %>%
    group_by(ssp, district_id, year) %>%
    summarise(
      log10_Pest = mean(log10_Pest, na.rm = TRUE),
      T          = mean(T,  na.rm = TRUE),
      RH         = mean(RH, na.rm = TRUE),
      n_active   = n(),
      .groups    = "drop"
    ) %>%
    filter(n_active >= 1)
}

# ---- (1) OLS hierarchical (yearly) ----
cat("\n========== OLS hierarchical (yearly, log10 P_est) ==========\n")
cat("SSP:", SSP_LABELS[TARGET_SSP], "(", TARGET_SSP, ")\n")

ols_panel <- load_pest_yearly(TARGET_SSP)
cat("OLS panel n =", nrow(ols_panel),
    "across", n_distinct(ols_panel$district_id), "districts\n")

m1_ols <- lm(log10_Pest ~ T + I(T^2),       data = ols_panel)
m2_ols <- lm(log10_Pest ~ T + I(T^2) + RH,  data = ols_panel)

ols_table <- tibble(
  Model    = c("M1: T + T^2", "M2: + RH"),
  R2       = c(summary(m1_ols)$r.squared, summary(m2_ols)$r.squared),
  Adj_R2   = c(summary(m1_ols)$adj.r.squared, summary(m2_ols)$adj.r.squared),
  delta_R2 = c(NA_real_, summary(m2_ols)$r.squared - summary(m1_ols)$r.squared),
  ANOVA_p  = c(NA_real_, anova(m1_ols, m2_ols)[["Pr(>F)"]][2])
)
print(ols_table)
write_csv(ols_table, file.path(TABLES_DIR, "ols_hierarchical.csv"))

b_T  <- coef(m2_ols)["T"]
b_T2 <- coef(m2_ols)["I(T^2)"]
T_opt_ols <- as.numeric(-b_T / (2 * b_T2))
cat(sprintf("OLS M2: beta_T = %.3f, beta_T2 = %.4f, T_opt = %.2f °C\n",
            b_T, b_T2, T_opt_ols))

# ---- (2) LMM (monthly) ----
cat("\n========== LMM (monthly panel, log10 P_est) ==========\n")

lmm_panel <- load_pest_panel(TARGET_SSP)
cat("LMM panel n =", nrow(lmm_panel),
    "across", n_distinct(lmm_panel$district_id), "districts\n")

# LMM1: random intercept
# LMM2: + random slope on T
lmm1 <- lmer(log10_Pest ~ T + I(T^2) + RH + (1 | district_id),
             data = lmm_panel, REML = FALSE)
lmm2 <- lmer(log10_Pest ~ T + I(T^2) + RH + (T | district_id),
             data = lmm_panel, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 2e5)))

lrt <- anova(lmm1, lmm2)

fit_metrics <- tibble(
  Model     = c("LMM1: T + T^2 + RH + (1|district)",
                "LMM2: + (T|district) random slope"),
  AIC       = c(AIC(lmm1), AIC(lmm2)),
  BIC       = c(BIC(lmm1), BIC(lmm2)),
  logLik    = c(as.numeric(logLik(lmm1)), as.numeric(logLik(lmm2))),
  delta_AIC = c(NA_real_, AIC(lmm1) - AIC(lmm2)),
  chi2      = c(NA_real_, lrt[["Chisq"]][2]),
  p_LRT     = c(NA_real_, lrt[["Pr(>Chisq)"]][2])
)
print(fit_metrics)
write_csv(fit_metrics, file.path(TABLES_DIR, "lmm_fit_metrics.csv"))

# REML for inference on fixed effects
lmm1_reml <- lmer(log10_Pest ~ T + I(T^2) + RH + (1 | district_id),
                  data = lmm_panel, REML = TRUE)
lmm2_reml <- lmer(log10_Pest ~ T + I(T^2) + RH + (T | district_id),
                  data = lmm_panel, REML = TRUE,
                  control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e5)))

fe_tab <- broom.mixed::tidy(lmm1_reml, effects = "fixed", conf.int = TRUE) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)
print(fe_tab)
write_csv(fe_tab, file.path(TABLES_DIR, "lmm_fixed_effects.csv"))

b_T_lmm  <- fixef(lmm1_reml)["T"]
b_T2_lmm <- fixef(lmm1_reml)["I(T^2)"]
T_opt_lmm <- as.numeric(-b_T_lmm / (2 * b_T2_lmm))
cat(sprintf("LMM1 (REML): beta_T = %.3f, beta_T2 = %.4f, T_opt = %.2f °C\n",
            b_T_lmm, b_T2_lmm, T_opt_lmm))

icc_val <- performance::icc(lmm1_reml)
cat("LMM1 ICC:\n"); print(icc_val)

# Nakagawa & Schielzeth (2013) marginal and conditional R^2
# Marginal R^2: variance explained by fixed effects only
# Conditional R^2: variance explained by fixed + random effects
r2_lmm1 <- performance::r2_nakagawa(lmm1_reml)
r2_lmm2 <- performance::r2_nakagawa(lmm2_reml)
cat("\nLMM1 (random intercept) Nakagawa R^2:\n")
print(r2_lmm1)
cat("\nLMM2 (random intercept + slope) Nakagawa R^2:\n")
print(r2_lmm2)

r2_table <- tibble(
  Model         = c("LMM1: random intercept",
                    "LMM2: random intercept + slope"),
  R2_marginal   = c(r2_lmm1$R2_marginal, r2_lmm2$R2_marginal),
  R2_conditional= c(r2_lmm1$R2_conditional, r2_lmm2$R2_conditional)
)
print(r2_table)
write_csv(r2_table, file.path(TABLES_DIR, "lmm_r2_nakagawa.csv"))

vc <- as.data.frame(VarCorr(lmm1_reml))
sigma_district <- sqrt(vc$vcov[vc$grp == "district_id" & vc$var1 == "(Intercept)"])
sigma_resid    <- sqrt(vc$vcov[vc$grp == "Residual"])
cat(sprintf("sigma_district = %.3f, sigma_resid = %.3f, ratio = %.2f\n",
            sigma_district, sigma_resid, sigma_district / sigma_resid))

# ---- (3) Cross-SSP T_opt ----
cat("\n========== Cross-SSP T_opt (LMM1 spec) ==========\n")

cross_ssp_results <- lapply(ALL_SSPS, function(s) {
  panel <- load_pest_panel(s)
  fit   <- lmer(log10_Pest ~ T + I(T^2) + RH + (1 | district_id),
                data = panel, REML = TRUE)
  bT  <- fixef(fit)["T"]
  bT2 <- fixef(fit)["I(T^2)"]
  r2  <- performance::r2_nakagawa(fit)
  tibble(
    SSP_code  = s,
    SSP_label = SSP_LABELS[s],
    n_obs     = nrow(panel),
    T_opt     = as.numeric(-bT / (2 * bT2)),
    ICC       = performance::icc(fit)$ICC_adjusted,
    R2_marg   = r2$R2_marginal,
    R2_cond   = r2$R2_conditional
  )
}) %>% bind_rows()
print(cross_ssp_results)
write_csv(cross_ssp_results, file.path(TABLES_DIR, "cross_ssp_topt.csv"))

# ---- (4) OLS vs LMM comparison ----
ols_lmm_comp <- tibble(
  Feature        = c("Bağımlı değişken", "Gözlem birimi", "n",
                     "İlçe etkisi",
                     "OLS R² (M2)",
                     "LMM R²_marginal (sabit etkiler)",
                     "LMM R²_conditional (sabit + rastgele)",
                     "ICC",
                     "T_optimum (°C)",
                     "Brière biyolojisiyle uyum (26-28°C)",
                     "Λ_import kaynaklı yapısal bağımlılık"),
  `OLS_yıllık`   = c("log10(P_est)", "İlçe-yıl", as.character(nrow(ols_panel)),
                     "Kontrol edilmemiş",
                     sprintf("%.3f", summary(m2_ols)$r.squared),
                     "—", "—", "—",
                     sprintf("%.1f", T_opt_ols),
                     "Yakın",
                     "Yok"),
  `LMM_aylık`    = c("log10(P_est)", "İlçe-yıl-ay", as.character(nrow(lmm_panel)),
                     "Random intercept",
                     "—",
                     sprintf("%.3f", r2_lmm1$R2_marginal),
                     sprintf("%.3f", r2_lmm1$R2_conditional),
                     sprintf("%.3f", performance::icc(lmm1_reml)$ICC_adjusted),
                     sprintf("%.1f", T_opt_lmm),
                     "Yakın",
                     "Yok")
)
print(ols_lmm_comp)
write_csv(ols_lmm_comp, file.path(TABLES_DIR, "ols_lmm_comparison.csv"))

cat("\n========== Done. Tables written to:", TABLES_DIR, "==========\n")