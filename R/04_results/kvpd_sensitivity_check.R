# =============================================================================
# kvpd_sensitivity_check.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Actually RUN the k_vpd sensitivity analysis that S2 File / the manuscript
#   claim was performed (k_vpd = 0.3, 0.5, 0.8 kPa^-1). Until now this was only
#   a to-do comment in the pipeline; the reported "38.3% / 38.3% / 38.8%" numbers
#   have no reproducible source. This script recomputes establishment probability
#   from the real district-month climate trajectories at each k_vpd value and
#   reports (i) whether the inter-district risk hierarchy is preserved and
#   (ii) the resulting maximum P_est and 50-year cumulative risk P_horizon.
#
# WHY RE-SIMULATE (rather than reuse the stored P_est)
#   k_vpd enters vector mortality:
#       mu_v(T,RH) = (1/lf(T)) * exp( k_vpd * (VPD(T,RH) - VPD_ref) )
#   which feeds lambda_i1 -> R0 -> P_est. The stored p_establishment_mean is
#   fixed at the baseline k_vpd = 0.5, so P_est must be recomputed from T and RH
#   for every alternative k_vpd. The monthly RDS carries temp_c and rh, so the
#   full mechanistic chain can be reconstructed exactly.
#
# INPUTS  (one monthly RDS per SSP scenario)
#   outputs/<ssp>/simulation/ctmc_spark_monthly_2025_2075_rep1000.rds
#     required columns: district_id, year, month, species, temp_c, rh,
#                       q_import_month   (importation side is k_vpd-independent)
#   The trait parameter table (species-specific Briere/quadratic coefficients,
#   VPD_ref, k_vpd) is read from the project data if available; otherwise the
#   canonical S1-File values are used (see TRAITS below).
#
# OUTPUTS
#   outputs/tables/tbl_kvpd_sensitivity.csv         (P_horizon & max P_est by k_vpd)
#   outputs/tables/tbl_kvpd_reference_check.csv     (RH=75% max-P_est check)
#   console summary + hierarchy-preservation verdict + baseline calibration
#
# MECHANISTIC FORMULAS  (verbatim from parameter_functions.R / ctmc_spark.R)
#   a(T)    = Briere(c,T0,Tm)
#   EIP(T)  = 1 / Briere_dev(c,T0,Tm)
#   lf(T)   = max( quadratic_unimodal(c,T0,Tm), lf_floor=0.25 )
#   VPD     = max( 0.6108*exp(17.27T/(T+237.3)) * (1 - RH/100), 0 )
#   mu_v    = max( (1/lf) * exp( k_vpd*(VPD - VPD_ref) ), 1e-6 )
#   lambda1 = m * a^2 * beta_vh * beta_hv * exp(-mu_v*EIP) / mu_v
#   R0      = lambda1 / gamma
#   P_est   = 1 - q(tau),  q from finite-threshold gambler's ruin (tau=30).
#             Because lambda_n = n*lambda1 and mu_n = n*gamma, the ratio
#             r_n = mu_n/lambda_n = 1/R0 is state-independent, so the general
#             birth-death extinction formula reduces to the closed form:
#                 q(tau) = (R0^-1 - R0^-tau)/(1 - R0^-tau),  R0 != 1
#                 q(tau) = 1 - 1/tau,                        R0  = 1
#             This exactly matches extinction_prob_bd() for the spark phase.
#
#   Jensen correction: individual EIP ~ LogNormal at fixed T with
#   sdlog(T)=exp(alpha+beta*(T-Tref)) calibrated to Chan & Johansson anchors
#   (25C:5-33d, 30C:2-15d), mean fixed by the Briere dev-rate. Inner MC over
#   n_mc EIP draws gives E[P_est]; set USE_JENSEN=TRUE to match the main outputs.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(here)
})

# ---- Configuration ----------------------------------------------------------
SSP_SCENARIOS <- c("ssp126", "ssp245", "ssp585")
SSP_LABELS    <- c(ssp126="SSP1-2.6", ssp245="SSP2-4.5", ssp585="SSP5-8.5")
K_VPD_VALUES  <- c(0.3, 0.5, 0.8)          # baseline = 0.5
MONTHLY_FILE  <- "ctmc_spark_monthly_2025_2075_rep1000.rds"

USE_JENSEN <- TRUE      # TRUE: nested-MC E[P_est] (matches main model). FALSE: deterministic.
N_MC       <- 2000      # inner EIP draws (only if USE_JENSEN)
SEED       <- 12345

# Fixed transmission parameters (S1 File)
BETA_VH <- 0.30; BETA_HV <- 0.33; GAMMA <- 0.20; M_RATIO <- 1.0
TAU <- 30L; LF_FLOOR <- 0.25; VPD_REF <- 1.0

DISTRICT_ORDER <- c("TUR.40.25_1","TUR.10.4_1","TUR.81.6_1","TUR.59.4_1","TUR.39.3_1")
DISTRICT_LABELS <- c("TUR.40.25_1"="Kartal","TUR.10.4_1"="Hopa",
                     "TUR.81.6_1"="Zonguldak","TUR.59.4_1"="Fethiye","TUR.39.3_1"="Egirdir")

# Canonical species trait coefficients (S1 File Table). Used if no project table is found.
# form: a_biting=Briere, eip_dev_rate=Briere, lifespan_temp=quadratic
TRAITS <- list(
  aegypti = list(
    a   = c(c = 2.71e-4, T0 = 14.67, Tm = 41.00),
    dev = c(c = 1.04e-4, T0 = 11.50, Tm = 38.97),
    lf  = c(c = 1.48e-1, T0 = 9.16,  Tm = 37.73)
  ),
  albopictus = list(
    a   = c(c = 1.93e-4, T0 = 10.25, Tm = 38.32),
    dev = c(c = 1.09e-4, T0 = 10.39, Tm = 43.05),
    lf  = c(c = 1.43e-1, T0 = 6.24,  Tm = 38.25)
  )
)

# ---- Thermal performance primitives (verbatim forms) ------------------------
briere <- function(T, c, T0, Tm) {
  out <- numeric(length(T)); ok <- is.finite(T) & T > T0 & T < Tm
  out[ok] <- c * T[ok] * (T[ok] - T0) * sqrt(Tm - T[ok]); pmax(out, 0)
}
quadratic_unimodal <- function(T, c, T0, Tm) {
  out <- numeric(length(T)); ok <- is.finite(T) & T > T0 & T < Tm
  out[ok] <- -c * (T[ok] - T0) * (T[ok] - Tm); pmax(out, 0)
}
vpd_kpa <- function(T, RH) {
  es <- 0.6108 * exp(17.27 * T / (T + 237.3)); pmax(es - es * RH / 100, 0)
}

# ---- EIP log-normal sdlog(T) calibration (Chan & Johansson anchors) ---------
.calibrate_sdlog <- function(dev_par, Tref = 27.5) {
  eip_mean <- function(T) { r <- briere(T, dev_par["c"], dev_par["T0"], dev_par["Tm"]); ifelse(r > 0, 1/r, Inf) }
  anchors <- list(list(T=25,q5=5,q95=33), list(T=30,q5=2,q95=15))
  z05 <- qnorm(0.05); z95 <- qnorm(0.95)
  qfun <- function(mean, s, p) if (!is.finite(mean) || mean <= 0) Inf else exp(log(mean) - 0.5*s^2 + s*qnorm(p))
  s0 <- function(q5,q95) log(q95/q5)/(z95 - z05)
  s25 <- s0(5,33); s30 <- s0(2,15)
  beta0 <- (log(s30)-log(s25))/(30-25); alpha0 <- log(s25) - beta0*(25 - Tref)
  obj <- function(th) { a<-th[1]; b<-th[2]; ss<-0
    for (an in anchors) { s<-exp(a+b*(an$T-Tref)); mE<-eip_mean(an$T)
      ss <- ss + (log(qfun(mE,s,0.05)/an$q5))^2 + (log(qfun(mE,s,0.95)/an$q95))^2 }
    ss }
  fit <- optim(c(alpha0,beta0), obj, method="BFGS", control=list(reltol=1e-10,maxit=200))
  list(alpha=fit$par[1], beta=fit$par[2], Tref=Tref)
}

# ---- Closed-form spark-phase establishment probability ----------------------
# R0 scalar -> P_est (finite-threshold gambler's ruin, tau=30)
p_est_from_R0 <- function(R0) {
  if (!is.finite(R0) || R0 <= 0) return(0)
  if (abs(R0 - 1) < 1e-12)       return(1 / TAU)
  q <- (R0^(-1) - R0^(-TAU)) / (1 - R0^(-TAU))
  min(max(1 - q, 0), 1)
}

# lambda_i1 for a scalar EIP value
.lambda1 <- function(aT, muv, eip) {
  surv <- if (is.finite(eip)) exp(-muv * eip) else 0
  M_RATIO * aT^2 * BETA_VH * BETA_HV * surv / muv
}

# Establishment probability for ONE district-month at a given k_vpd
p_est_one <- function(T, RH, species, k_vpd, use_jensen, n_mc, sdlog_fit) {
  tr <- TRAITS[[species]]
  aT  <- briere(T, tr$a["c"], tr$a["T0"], tr$a["Tm"])
  lf  <- max(quadratic_unimodal(T, tr$lf["c"], tr$lf["T0"], tr$lf["Tm"]), LF_FLOOR)
  muv <- max((1 / lf) * exp(k_vpd * (vpd_kpa(T, RH) - VPD_REF)), 1e-6)
  dev <- briere(T, tr$dev["c"], tr$dev["T0"], tr$dev["Tm"])
  eip_mean <- if (dev > 0) 1 / dev else Inf

  if (aT <= 0 || !is.finite(eip_mean)) return(0)   # outside thermal window

  if (!use_jensen) {
    return(p_est_from_R0(.lambda1(aT, muv, eip_mean) / GAMMA))
  }
  # Jensen-corrected: average P_est over log-normal EIP draws
  cf <- sdlog_fit[[species]]
  s  <- exp(cf$alpha + cf$beta * (T - cf$Tref))
  meanlog <- log(eip_mean) - 0.5 * s^2
  eip_draws <- rlnorm(n_mc, meanlog = meanlog, sdlog = s)
  mean(vapply(eip_draws, function(e) p_est_from_R0(.lambda1(aT, muv, e) / GAMMA), numeric(1)))
}

# ---- Path resolution --------------------------------------------------------
monthly_path <- function(ssp) {
  rel <- file.path("outputs", ssp, "simulation", MONTHLY_FILE)
  p <- tryCatch(here::here(rel), error = function(e) rel)
  if (!file.exists(p)) p <- rel
  p
}

# =============================================================================
# PART 1 — Full k_vpd sensitivity on REAL district-month trajectories
# =============================================================================
set.seed(SEED)
sdlog_fit <- lapply(TRAITS, function(tr) .calibrate_sdlog(tr$dev))

message("PART 1: recomputing P_est on real trajectories for k_vpd = ",
        paste(K_VPD_VALUES, collapse = ", "),
        if (USE_JENSEN) sprintf("  (Jensen-corrected, n_mc=%d)", N_MC) else "  (deterministic))")

sens_rows <- list(); calib_rows <- list()

for (ssp in SSP_SCENARIOS) {
  path <- monthly_path(ssp)
  if (!file.exists(path)) { warning(sprintf("[%s] missing: %s", ssp, path)); next }
  df <- readRDS(path)
  req <- c("district_id","species","temp_c","rh","q_import_month")
  if (!all(req %in% names(df))) stop(sprintf("[%s] missing columns: %s", ssp,
       paste(setdiff(req, names(df)), collapse=", ")))
  message(sprintf("  [%s] rows=%d", ssp, nrow(df)))

  for (k in K_VPD_VALUES) {
    pe <- vapply(seq_len(nrow(df)), function(i)
      p_est_one(df$temp_c[i], df$rh[i], df$species[i], k, USE_JENSEN, N_MC, sdlog_fit),
      numeric(1))
    pmm <- pmin(pmax(df$q_import_month * pe, 0), 1)   # p_month_major = q_import * P_est

    g <- tibble(district_id = df$district_id, pe = pe, pmm = pmm)
    agg <- g %>% group_by(district_id) %>%
      summarise(max_p_est = max(pe, na.rm = TRUE),
                P_horizon = 1 - prod(1 - pmm[!is.na(pmm)]), .groups = "drop") %>%
      mutate(ssp = ssp, ssp_label = unname(SSP_LABELS[ssp]), k_vpd = k,
             district = dplyr::recode(district_id, !!!DISTRICT_LABELS))
    sens_rows[[paste(ssp,k)]] <- agg

    # Baseline calibration: k=0.5 recomputed vs stored p_establishment_mean
    if (isTRUE(all.equal(k, 0.5)) && "p_establishment_mean" %in% names(df)) {
      act <- df$p_establishment_mean > 1e-12
      if (any(act)) calib_rows[[ssp]] <- tibble(
        ssp = ssp,
        corr = suppressWarnings(cor(pe[act], df$p_establishment_mean[act])),
        my_max = max(pe), stored_max = max(df$p_establishment_mean))
    }
  }
}

if (length(sens_rows) == 0) stop("No scenario files could be read.")

sens <- bind_rows(sens_rows) %>%
  mutate(ssp = factor(ssp, levels = SSP_SCENARIOS),
         district_id = factor(district_id, levels = DISTRICT_ORDER)) %>%
  arrange(ssp, district_id, k_vpd) %>%
  mutate(ssp = as.character(ssp), district_id = as.character(district_id))

# ---- Hierarchy-preservation verdict (per scenario, across k_vpd) ------------
hierarchy_ok <- sens %>%
  group_by(ssp, k_vpd) %>%
  arrange(desc(P_horizon), .by_group = TRUE) %>%
  summarise(order = paste(district, collapse = ">"), .groups = "drop") %>%
  group_by(ssp) %>%
  summarise(preserved = n_distinct(order) == 1,
            orders = paste(unique(order), collapse = "  |  "), .groups = "drop")

# ---- Write & print ----------------------------------------------------------
out_dir <- tryCatch(here::here("outputs","tables"), error=function(e) file.path("outputs","tables"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sens_out <- sens %>%
  transmute(Senaryo = ssp_label, Ilce = district,
            `k_vpd` = k_vpd, `maks P_est` = max_p_est, `P_horizon` = P_horizon)
readr::write_csv(sens_out, file.path(out_dir, "tbl_kvpd_sensitivity.csv"))

fmt <- function(x) formatC(x, format="e", digits=3)
cat("\n=========== k_vpd SENSITIVITY — P_horizon (real trajectories) ===========\n")
ph_wide <- sens %>% select(ssp_label, district, k_vpd, P_horizon) %>%
  tidyr::pivot_wider(names_from = k_vpd, values_from = P_horizon,
                     names_prefix = "k=")
print(as.data.frame(ph_wide %>% mutate(across(starts_with("k="), fmt))), row.names = FALSE)

cat("\n----------- max P_est by k_vpd -----------\n")
pe_wide <- sens %>% select(ssp_label, district, k_vpd, max_p_est) %>%
  tidyr::pivot_wider(names_from = k_vpd, values_from = max_p_est, names_prefix = "k=")
print(as.data.frame(pe_wide %>% mutate(across(starts_with("k="),
      ~formatC(.x, format="f", digits=4)))), row.names = FALSE)

cat("\n----------- HIERARCHY PRESERVATION -----------\n")
for (i in seq_len(nrow(hierarchy_ok))) {
  cat(sprintf("  %s: %s\n", SSP_LABELS[hierarchy_ok$ssp[i]],
      if (hierarchy_ok$preserved[i]) "PRESERVED across all k_vpd"
      else paste0("CHANGES -> ", hierarchy_ok$orders[i])))
}

if (length(calib_rows) > 0) {
  cat("\n----------- BASELINE CALIBRATION (k=0.5 recomputed vs stored) -----------\n")
  cb <- bind_rows(calib_rows)
  for (i in seq_len(nrow(cb)))
    cat(sprintf("  %s: corr=%.4f  my_max=%.4f  stored_max=%.4f\n",
        SSP_LABELS[cb$ssp[i]], cb$corr[i], cb$my_max[i], cb$stored_max[i]))
  cat("  (High corr confirms the recomputation matches the pipeline; residual\n",
      "   gap is the Jensen correction if USE_JENSEN differs from the stored run.)\n")
}

# =============================================================================
# PART 2 — Reference-condition check (reproduces the S2-style claim)
# =============================================================================
# S2 reports "at RH=75%, max P_est was 38.3/38.3/38.8% for k=0.3/0.5/0.8".
# This sweeps temperature at fixed RH=75% and reports max P_est per species,
# so the manuscript number can be checked against a reproducible value.
cat("\n=========== REFERENCE-CONDITION CHECK (RH = 75%, deterministic) ===========\n")
Tgrid <- seq(0, 45, by = 0.02)
ref_rows <- list()
for (sp in names(TRAITS)) {
  for (k in K_VPD_VALUES) {
    pe <- vapply(Tgrid, function(T) p_est_one(T, 75, sp, k, use_jensen = FALSE, N_MC, sdlog_fit), numeric(1))
    ref_rows[[paste(sp,k)]] <- tibble(species = sp, k_vpd = k,
                                      max_p_est = max(pe),
                                      T_at_max = Tgrid[which.max(pe)])
  }
}
ref <- bind_rows(ref_rows)
readr::write_csv(ref, file.path(out_dir, "tbl_kvpd_reference_check.csv"))
ref_wide <- ref %>% select(species, k_vpd, max_p_est) %>%
  tidyr::pivot_wider(names_from = k_vpd, values_from = max_p_est, names_prefix = "k=") %>%
  mutate(across(starts_with("k="), ~sprintf("%.1f%%", .x * 100)))
print(as.data.frame(ref_wide), row.names = FALSE)
cat("\nNOTE: compare these reproducible values against the number stated in S2 File.\n",
    "If they disagree, the S2 figure should be replaced with the value above\n",
    "(or the claim removed), since this is the reproducible reference computation.\n")

cat("\nDone. Tables written to:", out_dir, "\n")
