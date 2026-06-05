# ==========================================================
# R/03_models/parameter_functions.R
# Vector-trait parameter functions for dengue mechanistic models
# Focus: a(T), EIP(T), mu_v(T,RH) and derived lambda_local
#
# Key design choices (thesis-friendly):
# 1) Deterministic climate column detection (strict, prioritized)
# 2) Trait functions implemented in standard literature forms
#    - Brière / Quadratic thermal performance curves
#    - Mean EIP = 1 / development_rate(T)
#    - Humidity effect via VPD (vapour pressure deficit) multiplier
# 3) NO invented coefficients:
#    - coefficients must be supplied via a parameter table (CSV/RDS)
#
# Addendum (stochastic EIP):
# - Optional log-normal EIP draws with temperature-dependent sdlog(T)
#   calibrated to 5th–95th anchors:
#     * 25 °C: 5–33 days
#     * 30 °C: 2–15 days
#   (Anchors consistent with Chan & Johansson, 2012.)
# ==========================================================

## ----------------------------------------------------------
## 0) Dependencies (assume init.R already loads dplyr/purrr/etc.)
## ----------------------------------------------------------
# No library() calls here by design. Keep pure functions.

# Cache environment for EIP log-normal calibration (per R session)
.eip_logn_fit_env <- new.env(parent = emptyenv())

## ----------------------------------------------------------
## 1) Deterministic climate column detection
## ----------------------------------------------------------
.detect_first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

detect_temperature_col <- function(df) {
  candidates <- c("tas", "tas_mean", "t2m", "tmean", "temp_c", "temperature", "T", "temp")
  col <- .detect_first_present(names(df), candidates)
  if (is.na(col)) {
    stop(
      "Temperature column not found. Expected one of: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n  - ", paste(names(df), collapse = "\n  - "),
      call. = FALSE
    )
  }
  col
}

detect_rh_col <- function(df) {
  candidates <- c("hur", "hurs", "rh", "RH", "rel_humidity", "humidity")
  col <- .detect_first_present(names(df), candidates)
  if (is.na(col)) {
    stop(
      "Relative humidity (RH) column not found. Expected one of: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n  - ", paste(names(df), collapse = "\n  - "),
      call. = FALSE
    )
  }
  col
}

is_kelvin_temperature <- function(x) {
  stats::median(x, na.rm = TRUE) > 100
}

climate_long_to_wide <- function(df, key = c("district_id","year","month"),
                                 var_col = "variable", value_col = "value") {
  miss <- setdiff(c(key, var_col, value_col), names(df))
  if (length(miss) > 0) {
    stop("Climate long format is missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  df_wide <- df %>%
    dplyr::select(dplyr::all_of(key), dplyr::all_of(var_col), dplyr::all_of(value_col)) %>%
    dplyr::mutate(
      .var = as.character(.data[[var_col]]),
      .val = as.numeric(.data[[value_col]])
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key)), .data$.var) %>%
    dplyr::summarise(.val = mean(.data$.val, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = .data$.var, values_from = .data$.val)
  
  df_wide
}

standardize_climate_cols <- function(df, key = c("district_id","year","month")) {
  miss_key <- setdiff(key, names(df))
  if (length(miss_key) > 0) stop("Missing key columns: ", paste(miss_key, collapse = ", "), call. = FALSE)
  
  if (all(c("variable", "value") %in% names(df))) {
    df <- climate_long_to_wide(df, key = key, var_col = "variable", value_col = "value")
  }
  
  tcol <- detect_temperature_col(df)
  hcol <- detect_rh_col(df)
  
  Traw  <- as.numeric(df[[tcol]])
  RHraw <- as.numeric(df[[hcol]])
  
  if (is_kelvin_temperature(Traw)) {
    Traw <- Traw - 273.15
  }
  
  out <- df[, key, drop = FALSE]
  out$temp_c <- as.numeric(Traw)
  out$rh     <- as.numeric(RHraw)
  
  if (any(!is.finite(out$temp_c))) stop("temp_c has non-finite values after standardization.", call. = FALSE)
  if (any(is.na(out$rh)))          stop("rh has NA after standardization.", call. = FALSE)
  if (any(out$rh < 0 | out$rh > 100, na.rm = TRUE)) stop("rh must be in 0..100.", call. = FALSE)
  
  out
}

## ----------------------------------------------------------
## 2) Generic thermal performance curve forms
## ----------------------------------------------------------
briere <- function(T, c, T0, Tm) {
  out <- rep(0, length(T))
  ok  <- is.finite(T) & (T > T0) & (T < Tm)
  out[ok] <- c * T[ok] * (T[ok] - T0) * sqrt(Tm - T[ok])
  pmax(out, 0)
}

quadratic_unimodal <- function(T, c, T0, Tm) {
  out <- rep(0, length(T))
  ok  <- is.finite(T) & (T > T0) & (T < Tm)
  out[ok] <- -c * (T[ok] - T0) * (T[ok] - Tm)
  pmax(out, 0)
}

## ----------------------------------------------------------
## 3) Parameter table I/O and validation
## ----------------------------------------------------------
read_trait_params <- function(path) {
  if (!file.exists(path)) stop("Trait parameter file not found: ", path, call. = FALSE)
  ext <- tools::file_ext(path)
  
  if (tolower(ext) %in% c("rds","rda","rdata")) {
    tab <- readRDS(path)
  } else if (tolower(ext) %in% c("csv")) {
    tab <- readr::read_csv(path, show_col_types = FALSE)
  } else {
    stop("Unsupported parameter file extension: ", ext, call. = FALSE)
  }
  
  req  <- c("trait","form","c","T0","Tm")
  miss <- setdiff(req, names(tab))
  if (length(miss) > 0) {
    stop("Trait param table missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  tab
}

.get_trait_row <- function(tab, trait_name) {
  row <- tab[tab$trait == trait_name, , drop = FALSE]
  if (nrow(row) != 1) {
    stop(
      "Trait '", trait_name, "' must have exactly 1 row in parameter table. Found: ", nrow(row),
      call. = FALSE
    )
  }
  row
}

## ----------------------------------------------------------
## 4) Trait functions: a(T), mean EIP(T), mu_v(T,RH)
## ----------------------------------------------------------
a_of_T <- function(T, tab) {
  r    <- .get_trait_row(tab, "a_biting")
  form <- tolower(r$form[[1]])
  
  if (form == "briere") {
    briere(T, c = r$c[[1]], T0 = r$T0[[1]], Tm = r$Tm[[1]])
  } else if (form == "quadratic") {
    quadratic_unimodal(T, c = r$c[[1]], T0 = r$T0[[1]], Tm = r$Tm[[1]])
  } else {
    stop("Unsupported form for a_biting: ", r$form[[1]], call. = FALSE)
  }
}

# Mean EIP(T) in days: EIP = 1 / dev_rate(T)
eip_mean_of_T <- function(T, tab) {
  r    <- .get_trait_row(tab, "eip_dev_rate")
  form <- tolower(r$form[[1]])
  
  rate <- if (form == "briere") {
    briere(T, c = r$c[[1]], T0 = r$T0[[1]], Tm = r$Tm[[1]])
  } else if (form == "quadratic") {
    quadratic_unimodal(T, c = r$c[[1]], T0 = r$T0[[1]], Tm = r$Tm[[1]])
  } else {
    stop("Unsupported form for eip_dev_rate: ", r$form[[1]], call. = FALSE)
  }
  
  ifelse(rate > 0, 1 / rate, Inf)
}

vpd_kpa_from_TRH <- function(T, RH) {
  es  <- 0.6108 * exp(17.27 * T / (T + 237.3))
  ea  <- es * (RH / 100)
  vpd <- es - ea
  pmax(vpd, 0)
}

mu_v_of_TRH <- function(T, RH, tab,
                        lf_floor = 0.25) {   # gün; duyarlılık analizinde test et
  # ------------------------------------------------------------------
  # Trait tablosundaki "lifespan_temp" satırı, yetişkin YAŞAM SÜRESİ
  # lf(T) için fit edilmiş parametreleri içerir (birim: gün).
  # Mortalite hızı mu_v = 1 / lf(T) olarak türetilir.
  #
  # Biyolojik beklenti:
  #   lf(T) : tek tepeli (quadratic), optimum ~27-28 °C (aegypti)
  #   mu_v(T) = 1/lf(T) : optimumda MİNİMUM, uçlarda yüksek  <- U-şeklinde
  #
  # Numerik stabilizasyon:
  #   T0/Tm dışında lf(T) sıfıra yaklaşır; bu çok kısa yaşam süresi ve
  #   çok yüksek mortalite anlamına gelir. NA / Inf üretmemek için
  #   lf'ye alt sınır (lf_floor, varsayılan 0.25 gün) uygulanır.
  #   Bu, modelin termal sınırlar dışında da çökmeden devam etmesini sağlar.
  # ------------------------------------------------------------------
  
  rT    <- .get_trait_row(tab, "lifespan_temp")
  formT <- tolower(rT$form[[1]])
  
  lfT <- if (formT == "quadratic") {
    quadratic_unimodal(T, c = rT$c[[1]], T0 = rT$T0[[1]], Tm = rT$Tm[[1]])
  } else if (formT == "briere") {
    briere(T, c = rT$c[[1]], T0 = rT$T0[[1]], Tm = rT$Tm[[1]])
  } else {
    stop("Unsupported form for lifespan_temp: ", rT$form[[1]], call. = FALSE)
  }
  
  # Alt sınır: termal pencere dışında çok yüksek mortalite, ama Inf/NA değil
  lfT <- pmax(lfT, lf_floor)
  muT <- 1 / lfT   # 1/day
  
  rvpd <- .get_trait_row(tab, "vpd_ref")
  rk   <- .get_trait_row(tab, "k_vpd")
  
  if (tolower(rvpd$form[[1]]) != "scalar" || tolower(rk$form[[1]]) != "scalar") {
    stop("vpd_ref and k_vpd must have form='scalar' in parameter table.", call. = FALSE)
  }
  
  # VPD nem düzeltme çarpanı: düşük nem -> yüksek VPD -> yüksek mortalite
  vpd  <- vpd_kpa_from_TRH(T, RH)
  mult <- exp(rk$c[[1]] * (vpd - rvpd$c[[1]]))
  
  mu <- muT * mult
  
  # mu_floor: çok nadir ama VPD çarpanı sıfıra yaklaşırsa stabilite için
  mu <- pmax(mu, 1e-6)
  
  if (any(!is.finite(mu))) stop("mu_v produced non-finite values. Check parameter table and inputs.", call. = FALSE)
  if (any(mu < 0, na.rm = TRUE)) stop("mu_v produced negative values.", call. = FALSE)
  
  mu
}

## ----------------------------------------------------------
## 4.5) NEW: Log-normal EIP helpers (sdlog(T) calibration + caching)
## ----------------------------------------------------------

# Quantile of LogNormal given mean (E[X]) and sdlog
# meanlog is implied by mean: mean = exp(meanlog + 0.5*sdlog^2)
.lognorm_quantile_from_mean <- function(mean, sdlog, p) {
  if (!is.finite(mean) || mean <= 0) return(Inf)
  z <- stats::qnorm(p)
  exp(log(mean) - 0.5 * sdlog^2 + sdlog * z)
}

# Calibrate sdlog(T) with a simple parametric form:
#   sdlog(T) = exp(alpha + beta*(T - Tref))
# against thesis anchors (5th–95th ranges at 25°C and 30°C).
# Mean EIP(T) is fixed by eip_mean_of_T(T, tab) (NOT refit).
.calibrate_eip_logn_sdlog <- function(tab,
                                      anchors = list(
                                        list(T = 25, q5 = 5,  q95 = 33),
                                        list(T = 30, q5 = 2,  q95 = 15)
                                      ),
                                      Tref = 27.5) {
  
  z05 <- stats::qnorm(0.05)
  z95 <- stats::qnorm(0.95)
  
  # Initial guess for sdlog from ratio q95/q5 (independent of mean)
  s_from_ratio <- function(q5, q95) log(q95 / q5) / (z95 - z05)
  s25_init <- s_from_ratio(anchors[[1]]$q5, anchors[[1]]$q95)
  s30_init <- s_from_ratio(anchors[[2]]$q5, anchors[[2]]$q95)
  
  T1 <- anchors[[1]]$T; T2 <- anchors[[2]]$T
  beta0  <- (log(s30_init) - log(s25_init)) / (T2 - T1)
  alpha0 <- log(s25_init) - beta0 * (T1 - Tref)
  
  obj <- function(theta) {
    alpha <- theta[1]; beta <- theta[2]
    ss <- 0
    for (a in anchors) {
      T   <- a$T
      sT  <- exp(alpha + beta * (T - Tref))
      mE  <- eip_mean_of_T(T, tab)
      Q5  <- .lognorm_quantile_from_mean(mE, sT, 0.05)
      Q95 <- .lognorm_quantile_from_mean(mE, sT, 0.95)
      e5  <- log(Q5 / a$q5)
      e95 <- log(Q95 / a$q95)
      ss  <- ss + e5*e5 + e95*e95
    }
    ss
  }
  
  fit <- stats::optim(
    par     = c(alpha0, beta0),
    fn      = obj,
    method  = "BFGS",
    control = list(reltol = 1e-10, maxit = 200)
  )
  
  list(
    alpha = fit$par[1],
    beta  = fit$par[2],
    Tref  = Tref,
    value = fit$value,
    convergence = fit$convergence
  )
}

# Public: compute (meanlog, sdlog) for EIP LogNormal at temperature T
# Cache the calibration coefficients per session to avoid repeated optim().
get_eip_logn_params <- function(T, tab, force_recalibrate = FALSE) {
  if (isTRUE(force_recalibrate) || is.null(.eip_logn_fit_env$coef)) {
    .eip_logn_fit_env$coef <- .calibrate_eip_logn_sdlog(tab)
  }
  cf <- .eip_logn_fit_env$coef
  
  sdlog <- exp(cf$alpha + cf$beta * (T - cf$Tref))
  mean_eip <- eip_mean_of_T(T, tab)
  
  # Outside valid thermal range -> mean_eip may be Inf; treat as "no progression"
  if (!is.finite(mean_eip) || mean_eip <= 0) {
    return(list(meanlog = NA_real_, sdlog = as.numeric(sdlog)))
  }
  
  meanlog <- log(mean_eip) - 0.5 * sdlog^2
  list(meanlog = as.numeric(meanlog), sdlog = as.numeric(sdlog))
}

## ----------------------------------------------------------
## 5) Derived local birth intensity for CTMC
##     CORRECTED Ross-Macdonald formula:
##     lambda_i(t) = i * (m * a(T)^2 * beta_vh * beta_hv * exp(-mu_v * EIP)) / mu_v
##     Optional: EIP ~ LogNormal(meanlog(T), sdlog(T)) via Monte Carlo
## ----------------------------------------------------------
make_lambda_local <- function(n_infected, T, RH, beta_vh, beta_hv, m, tab_traits,
                              use_stochastic_EIP = FALSE,
                              n_mc = 2000L,
                              seed = NULL,
                              return_draws = FALSE) {
  
  if (length(T) != 1 || length(RH) != 1) {
    stop("make_lambda_local expects scalar T and RH for a given month (replicate upstream).", call. = FALSE)
  }
  
  n_infected <- as.integer(n_infected)
  if (any(n_infected < 0)) stop("n_infected must be >=0.", call. = FALSE)
  
  aT  <- a_of_T(T, tab_traits)            # 1/day
  muv <- mu_v_of_TRH(T, RH, tab_traits)   # 1/day
  
  # ============================================================
  # Ross-Macdonald formula:
  # λ = n_infected × (m·a²·β_vh·β_hv·S) / μ_v
  # 
  # Components:
  #   - m: mosquito/human ratio (vector abundance)
  #   - a²: biting rate squared (mosquito bites twice: 
  #         once to acquire virus, once to transmit)
  #   - β_vh: vector-to-human transmission probability
  #   - β_hv: human-to-vector transmission probability
  #   - S: survival through extrinsic incubation period
  #   - μ_v: vector mortality rate (population equilibrium)
  #
  # When use_stochastic_EIP = TRUE and return_draws = TRUE:
  #   Returns a VECTOR of n_mc λ values (one per individual EIP draw).
  #   The caller (compute_ctmc_month_spark) computes P_est for each draw
  #   and averages: E[P_est(EIP_i)].  This is the Jensen-correct approach
  #   described in Appendix 1, Section 3 (Eq. Ek1.12).
  #
  # When return_draws = FALSE (default, backward-compatible):
  #   Returns a scalar mean λ (old behaviour).
  # ============================================================
  
  if (isTRUE(use_stochastic_EIP)) {
    pars <- get_eip_logn_params(T, tab_traits, force_recalibrate = FALSE)
    
    if (!is.finite(pars$meanlog)) {
      # Outside thermal range: no EIP progression -> zero transmission
      if (isTRUE(return_draws)) {
        return(rep(0, as.integer(n_mc)))
      } else {
        return(0)
      }
    }
    
    if (!is.null(seed)) set.seed(seed)
    n_mc <- as.integer(n_mc)
    if (!is.finite(n_mc) || n_mc < 100L) {
      stop("n_mc must be an integer >= 100 for stochastic EIP.", call. = FALSE)
    }
    
    eip_draws <- stats::rlnorm(n_mc, meanlog = pars$meanlog, sdlog = pars$sdlog)
    surv_draws <- exp(-muv * eip_draws)
    
    # Per-draw λ values (one for each individual EIP realisation)
    lam_draws <- n_infected * (m * aT^2 * beta_vh * beta_hv * surv_draws) / muv
    lam_draws[!is.finite(lam_draws)] <- 0
    if (any(lam_draws < 0, na.rm = TRUE)) {
      stop("lambda_local negative: check parameterization.", call. = FALSE)
    }
    
    if (isTRUE(return_draws)) {
      return(lam_draws)
    } else {
      # Backward-compatible: return scalar mean
      return(mean(lam_draws))
    }
    
  } else {
    # Deterministic EIP: single scalar λ
    eip_mean <- eip_mean_of_T(T, tab_traits)
    surv <- if (!is.finite(eip_mean)) 0.0 else exp(-muv * eip_mean)
    
    lam <- n_infected * (m * aT^2 * beta_vh * beta_hv * surv) / muv
    lam[!is.finite(lam)] <- 0
    if (any(lam < 0, na.rm = TRUE)) {
      stop("lambda_local negative: check parameterization.", call. = FALSE)
    }
    
    if (isTRUE(return_draws)) {
      # Deterministic: return single-element vector for uniform interface
      return(lam)
    }
    return(lam)
  }
}

## ----------------------------------------------------------
## 6) Convenience: build a validated trait table template
## ----------------------------------------------------------
write_trait_param_template <- function(path = "data_processed/trait_params_template.csv") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  template <- data.frame(
    trait = c("a_biting", "eip_dev_rate", "lifespan_temp", "vpd_ref", "k_vpd"),
    form  = c("briere", "briere", "quadratic", "scalar", "scalar"),
    c     = c(NA, NA, NA, NA, NA),
    T0    = c(NA, NA, NA, NA, NA),
    Tm    = c(NA, NA, NA, NA, NA),
    unit  = c("1/day", "1/day", "day", "kPa", "per_kPa"),
    notes = c(
      "Biting rate a(T); c unit: 1/(day*degC^3.5); fill c,T0,Tm from your selected source",
      "EIP development rate; mean EIP=1/rate; fill c,T0,Tm from your selected source",
      "Adult lifespan lf(T) in DAYS (quadratic/briere); mu_v = 1/lf computed in mu_v_of_TRH(); fill c,T0,Tm from Mordecai 2017",
      "Reference VPD for humidity multiplier (kPa)",
      "Sensitivity of mortality to VPD difference"
    ),
    source = c(
      "Trait table source (e.g., fitted thermal curves)",
      "Trait table source (e.g., fitted thermal curves)",
      "Trait table source (e.g., fitted thermal curves / hazard model)",
      "Calibration target / literature",
      "Calibration target / literature"
    ),
    stringsAsFactors = FALSE
  )
  
  readr::write_csv(template, path, na = "")
  invisible(path)
}


message("✓ parameter_functions.R loaded (Ross-Macdonald formula + return_draws for Jensen-correct EIP MC)")
message("  make_lambda_local args: ", paste(names(formals(make_lambda_local)), collapse = ", "))