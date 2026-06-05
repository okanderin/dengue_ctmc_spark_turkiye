# ==========================================================
# R/03_models/ctmc_spark.R
# Spark-phase CTMC (birth-death) + species-aware traits
# - Climate: long format (variable/value) -> temp_c, rh
# - Importation: monthly expected imports (lambda_import)
# - Species rule: Artvin + Zonguldak -> aegypti; others -> albopictus
# - Outputs: monthly + yearly + horizon risk
# - Flexible "major" definitions via major_rule
# - SSP-AWARE: reads from DIR_PROCESSED_SSP, writes to DIR_OUTPUT_SSP
#
# CORRECTED VERSION: beta_hv parameter added to lambda calculation
# ==========================================================

## ----------------------------------------------------------
## 0) Initialize project environment

### 1) Kurulum ve İnputların Hazırlanması

#### init.R ve parameter_functions.R içe aktarılıyor. İkincisi; trait parametrelerini okuyacak (read_trait_params) ve ay–ilçe bazında yerel bulaş hızını (i=1 için) üretecek make_lambda_local(...) fonksiyonunu sağlamalı.

#### fp_import yolu yoksa işlenmiş aylık ithal baskısı dosyası (importation_pressure_monthly_2025_2075.rds) varsayılan tahmin edilen yoldan bulunuyor.


#### Neden? Modelin iki ana girdisi var: iklim (sıcaklık, RH) ve ithal vaka baskısı. Trait parametreleri tür bazında ısırma, mortalite, EIP gibi iklim-duyarlı bileşenleri sağlar; bulaş hızını doğrudan etkiler.
## ----------------------------------------------------------

# =============================================================================
# CTMC spark-phase model for dengue establishment in Turkey
# This script implements the early stochastic "spark" phase described in:
#   - Section 3.4 (Stochastic model, CTMC birth–death process)
#   - Section 3.5.3 (Monthly CTMC solution)
# =============================================================================

# Initialize project environment:
# - load directory paths (e.g., DIR_PROCESSED)
# - attach required packages
# - set global options used across the dengue CTMC spark-phase pipeline
source("R/01_setup/init.R")

# Load biological/trait-based parameter functions:
# - read_trait_params(): temperature- and humidity-dependent parameter sets
#   for a(T), μ_v(T,RH), EIP(T) etc. (Section 3.5.2.1; Eq. (1))
# - make_lambda_local(): construct the local transmission rate λ_local(I,t)
#   combining trait functions into the CTMC birth rate (Eq. (1)-(2))
source("R/03_models/parameter_functions.R")

# --- Importation pressure input file ----------------------------------------
# The spark-phase CTMC model requires pre-computed monthly importation pressure:
#   - Λ_import,month and q_import,month (Eq. (6)-(7); Section 3.4.4)
#   - indexed by district-month-year (Section 3.5.1)
# If fp_import was not defined earlier (e.g. in R/01_setup/paths.R),
# attempt to use the default processed file:
if (!exists("fp_import")) {
  # SSP-aware: importation data is now under DIR_PROCESSED_SSP
  fp_import_guess <- file.path(DIR_PROCESSED_SSP, "importation_pressure_monthly_2025_2075.rds")
  if (!file.exists(fp_import_guess)) {
    stop(
      "fp_import is not defined and the guessed file does not exist: ",
      fp_import_guess,
      "\nFix: define fp_import in R/01_setup/paths.R (preferred) ",
      "or ensure SSP_SCENARIO env var is set and data pipeline has run."
    )
  }
  fp_import <- fp_import_guess
}

# Null-coalescing helper:
# a %||% b returns 'a' if it is not NULL, otherwise returns 'b'.
# Used to provide sensible defaults for optional arguments (e.g., parameter sets).
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- Definition of "major outbreak" for the spark-phase BD process -----------------
# major_rule$type ∈ {"establishment", "hit_K", "sustain_days", "secondary_ge_X"}
#
# - "establishment" (default): major outbreak ≡ single imported infection
#   leads to non-extinction of the BD process; probability P_est = 1 − q
#   (Section 3.4.2–3.4.3; Eq. (3)).
# - "hit_K": BD process reaches I ≥ K at least once before extinction;
#   an analytically tractable threshold definition (hitting probability).
# - "sustain_days": BD process survives ≥ D days without extinction
#   (estimated via Gillespie simulations).
# - "secondary_ge_X": total number of local secondary cases ≥ X
#   (estimated via Gillespie simulations).
#
# Rationale: The literature uses different operational thresholds for
# "major outbreak" depending on context. This structure allows alternative
# definitions to be implemented on the same BD/CTMC backbone.



.major_rule_validate <- function(major_rule) {
  if (is.null(major_rule) || is.null(major_rule$type)) {
    return(list(type = "establishment"))  # backward-compatible default
  }
  t <- match.arg(major_rule$type, c("establishment","hit_K","sustain_days","secondary_ge_X"))
  out <- list(type = t)
  if (t == "hit_K") {
    K <- as.integer(major_rule$K %||% 20L)
    if (is.na(K) || K < 2L) stop("major_rule$K must be ≥ 2 for hit_K.")
    out$K <- K
  } else if (t == "sustain_days") {
    D <- as.numeric(major_rule$D_days %||% 14)
    n_path <- as.integer(major_rule$n_path %||% 2000L)
    if (!is.finite(D) || D <= 0) stop("major_rule$D_days must be > 0.")
    if (n_path < 100L) warning("n_path is small; increase to >1000 for stability.")
    out$D_days <- D; out$n_path <- n_path; out$seed <- major_rule$seed %||% NULL
  } else if (t == "secondary_ge_X") {
    X <- as.integer(major_rule$X %||% 20L)
    n_path <- as.integer(major_rule$n_path %||% 2000L)
    if (is.na(X) || X < 1L) stop("major_rule$X must be ≥ 1.")
    out$X <- X; out$n_path <- n_path; out$seed <- major_rule$seed %||% NULL
  }
  out
}

# --- Gillespie BD simülasyonu: I(0)=1, import KAPALI (yalnızca yerel) ---
.sim_bd_once <- function(lambda_local_i1, gamma_day, D_target = NULL, X_target = NULL,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- 1L               # I(0)=1
  t <- 0.0
  sec_cases <- 0L       # toplam yerel ikincil vaka sayısı
  while (n > 0L) {
    lam <- n * lambda_local_i1
    mu  <- n * gamma_day
    tot <- lam + mu
    if (tot <= 0) break
    dt <- rexp(1L, rate = tot)
    t  <- t + dt
    # Hedef süre sağlandıysa erken bitir
    if (!is.null(D_target) && t >= D_target) return(list(hit = TRUE, t = t, sec = sec_cases))
    # Olay seçimi
    if (runif(1) < lam / tot) {
      n <- n + 1L
      sec_cases <- sec_cases + 1L
      if (!is.null(X_target) && sec_cases >= X_target) return(list(hit = TRUE, t = t, sec = sec_cases))
    } else {
      n <- n - 1L
    }
  }
  list(hit = FALSE, t = t, sec = sec_cases)
}

# Monte Carlo estimator of P_major under duration or secondary-case thresholds.
.sim_bd_prob <- function(lambda_local_i1, gamma_day,
                         D_target = NULL, X_target = NULL,
                         n_path = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  hits <- 0L
  for (r in seq_len(n_path)) {
    s_here <- if (is.null(seed)) NULL else as.integer((seed + r) %% .Machine$integer.max)
    res <- .sim_bd_once(lambda_local_i1, gamma_day, D_target, X_target, seed = s_here)
    if (isTRUE(res$hit)) hits <- hits + 1L
  }
  hits / n_path
}

## ----------------------------------------------------------
## 1) Utilities: IO + asserts
## ----------------------------------------------------------
read_rds_checked <- function(path) {
  if (is.null(path) || length(path) != 1L || !is.character(path)) {
    stop("Invalid path argument.", call. = FALSE)
  }
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  readRDS(path)
}

assert_has_cols <- function(df, req, name = deparse(substitute(df))) {
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) {
    stop("Missing required columns in ", name, ": ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_unique_key <- function(df, key, name = deparse(substitute(df))) {
  assert_has_cols(df, key, name)
  n_all  <- nrow(df)
  n_uniq <- dplyr::n_distinct(df[, key, drop = FALSE])
  if (n_all != n_uniq) {
    stop("Key is not unique in ", name, ". Expected unique rows by: ",
         paste(key, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_numeric_nonneg <- function(x, nm) {
  if (!is.numeric(x)) stop(nm, " must be numeric.", call. = FALSE)
  if (any(!is.finite(x))) stop(nm, " has NA/non-finite.", call. = FALSE)
  if (any(x < 0)) stop(nm, " has negative values.", call. = FALSE)
  invisible(TRUE)
}

## ----------------------------------------------------------
## 2) Birth–death extinction probability (absorption at 0 vs tau)
## ----------------------------------------------------------

## ----------------------------------------------------------
## Birth–death extinction probability (absorption at 0 vs tau)
## q_extinction(i0) = P( I_t süreci 0'a ulaşsın, tau'ya ulaşmadan )
## ----------------------------------------------------------
extinction_prob_bd <- function(lambda_n, mu_n, tau, i0 = 1L) {
  # lambda_n, mu_n: vektörler, n = 1..(tau-1) için doğum/ölüm hızları
  # tau: üst sınır (absorbing) durum
  # i0: başlangıç enfekte sayısı (1..tau-1)
  
  if (length(lambda_n) != (tau - 1) || length(mu_n) != (tau - 1)) {
    stop("lambda_n and mu_n must have length tau-1.", call. = FALSE)
  }
  if (i0 < 1 || i0 > (tau - 1)) {
    stop("i0 must be in 1..(tau-1).", call. = FALSE)
  }
  
  # Eğer i0'dan itibaren herhangi bir durumda doğum hızı 0 veya altındaysa,
  # süreç yukarı ilerleyemeyeceği için 0'a sönme kaçınılmazdır.
  if (any(lambda_n[i0:(tau - 1)] <= 0)) return(1)
  
  # r_n = mu_n / lambda_n oranlarının log'ları
  log_r <- log(mu_n) - log(lambda_n)
  # B_k = prod_{j=1}^k r_j => log(B_k) = cumsum(log_r)
  log_B <- cumsum(log_r)
  
  logsumexp <- function(v) {
    m <- max(v)
    m + log(sum(exp(v - m)))
  }
  
  # Toplam payda: S_0 + S_1 + ... + S_{tau-1},
  # burada S_0 = 1, S_k = B_k (k>=1)
  # log( S_0 + ... + S_{tau-1} ) = logsumexp( c(log(1), log_B[1:(tau-1)]) )
  log_sum_all <- logsumexp(c(0, log_B[1:(tau - 1)]))
  
  # 0..(i0-1) arası prefix toplamı: S_0 + ... + S_{i0-1}
  if (i0 == 1L) {
    # S_0 = 1
    log_sum_prefix <- 0 # log(1)
  } else {
    log_sum_prefix <- logsumexp(c(0, log_B[1:(i0 - 1)]))
  }
  
  # q_i0 = 1 - ( S_0 + ... + S_{i0-1} ) / ( S_0 + ... + S_{tau-1} )
  q <- 1 - exp(log_sum_prefix - log_sum_all)
  
  # Sayısal güvenlik: [0,1] aralığına projeksiyon
  q <- min(max(q, 0), 1)
  q
}





## ----------------------------------------------------------
## 3) Load importation (monthly) -> lambda_import (expected imports/month)
## ----------------------------------------------------------
load_importation_monthly <- function(path, scenario_use = c("main","low","high")) {
  scenario_use <- match.arg(scenario_use)
  df <- read_rds_checked(path)
  
  # --- Choose ONE pi_scenario
  if ("pi_scenario" %in% names(df)) {
    df <- df %>%
      dplyr::filter(.data$pi_scenario == scenario_use) %>%
      dplyr::select(-.data$pi_scenario)
  }
  
  key <- c("district_id", "year", "month")
  assert_has_cols(df, key, "df_import")
  
  if (any(df$month < 1 | df$month > 12, na.rm = TRUE)) {
    stop("df_import has invalid month values (must be 1..12).", call. = FALSE)
  }
  
  # Construct lambda_import (monthly expected imports)
  if ("lambda_import" %in% names(df)) {
    df$lambda_import <- as.numeric(df$lambda_import)
  } else if ("expected_imported_cases_per_month" %in% names(df)) {
    df$lambda_import <- as.numeric(df$expected_imported_cases_per_month)
  } else if (all(c("lambda_import_per_day", "days_in_month") %in% names(df))) {
    df$lambda_import <- as.numeric(df$lambda_import_per_day) * as.numeric(df$days_in_month)
  } else {
    stop(
      "Cannot construct monthly lambda_import. Need one of:\n",
      "  - lambda_import\n",
      "  - expected_imported_cases_per_month\n",
      "  - lambda_import_per_day + days_in_month\n\n",
      "Available columns:\n  - ",
      paste(names(df), collapse = "\n  - "),
      call. = FALSE
    )
  }
  assert_numeric_nonneg(df$lambda_import, "df_import$lambda_import")
  
  # ---- Prefer q_import_month when available (Bernoulli importation)
  if ("q_import_month" %in% names(df)) {
    df$q_import_month <- as.numeric(df$q_import_month)
  } else if (all(c("lambda_import_per_day", "d_days", "days_in_month") %in% names(df))) {
    # If dataset provides a window length d_days and daily intensity:
    df$q_import_month <- 1 - exp(-as.numeric(df$lambda_import_per_day) * pmin(as.numeric(df$d_days), as.numeric(df$days_in_month)))
  } else if ("lambda_import" %in% names(df) && "days_in_month" %in% names(df)) {
    # fallback: interpret lambda_import as expected/month and approximate daily by /days (3-day window)
    df$q_import_month <- 1 - exp(-(as.numeric(df$lambda_import) / as.numeric(df$days_in_month)) * pmin(3, as.numeric(df$days_in_month)))
  } else {
    # as a last resort, set to 0
    df$q_import_month <- 0
  }
  assert_numeric_nonneg(df$q_import_month, "df_import$q_import_month")
  df$q_import_month <- pmin(pmax(df$q_import_month, 0), 1)
  
  # days_in_month: if missing or NA, compute
  if (!("days_in_month" %in% names(df))) {
    df$days_in_month <- as.integer(lubridate::days_in_month(lubridate::make_date(df$year, df$month, 1)))
  } else {
    df$days_in_month <- as.integer(df$days_in_month)
    idx_na <- which(is.na(df$days_in_month))
    if (length(idx_na) > 0) {
      df$days_in_month[idx_na] <- as.integer(lubridate::days_in_month(
        lubridate::make_date(df$year[idx_na], df$month[idx_na], 1)
      ))
    }
  }
  
  has_prov <- "province_name" %in% names(df)
  has_dist <- "district_name" %in% names(df)
  
  n_all  <- nrow(df)
  n_uniq <- dplyr::n_distinct(df[, key, drop = FALSE])
  
  if (n_all != n_uniq) {
    message("ℹ df_import has non-unique keys AFTER scenario filter. Aggregating by SUM(lambda_import) and OR-combining q_import_month.")
    df <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
      dplyr::summarise(
        lambda_import  = sum(.data$lambda_import, na.rm = TRUE),
        # 'at least one introduction' probability combined by OR:
        q_import_month = 1 - prod(1 - pmin(pmax(.data$q_import_month, 0), 1), na.rm = TRUE),
        days_in_month  = dplyr::first(.data$days_in_month),
        province_name  = if (has_prov) dplyr::first(stats::na.omit(as.character(.data$province_name))) else NA_character_,
        district_name  = if (has_dist) dplyr::first(stats::na.omit(as.character(.data$district_name))) else NA_character_,
        .groups = "drop"
      )
  } else {
    keep <- c(key, "lambda_import", "q_import_month", "days_in_month",
              intersect(names(df), c("province_name", "district_name")))
    df <- df %>% dplyr::select(dplyr::all_of(keep))
  }
  
  assert_unique_key(df, key, "df_import")
  df
}

## ----------------------------------------------------------
## 4) Load climate (monthly long) -> temp_c + rh
##    Expected columns: district_id, year, month, variable, value
##    Deterministic mapping:
##      temperature variable: tas (preferred) else temp_c/tmean/temperature
##      RH variable: hur (preferred) else rh/humidity/rel_humidity
## ----------------------------------------------------------
load_climate_monthly_long <- function(path) {
  df <- read_rds_checked(path)
  
  req <- c("district_id", "year", "month", "variable", "value")
  assert_has_cols(df, req, "df_climate_raw")
  
  if (any(df$month < 1 | df$month > 12, na.rm = TRUE)) {
    stop("df_climate_raw has invalid month values (must be 1..12).", call. = FALSE)
  }
  
  df$variable <- as.character(df$variable)
  df$value <- as.numeric(df$value)
  
  temp_vars <- c("tas", "temp_c", "tmean", "temperature")
  rh_vars   <- c("hur", "rh", "humidity", "rel_humidity")
  
  df2 <- df %>% dplyr::filter(.data$variable %in% c(temp_vars, rh_vars))
  if (nrow(df2) == 0) {
    stop(
      "Climate long file has no recognized temperature/RH variables.\n",
      "Expected temp one of: ", paste(temp_vars, collapse = ", "), "\n",
      "Expected RH   one of: ", paste(rh_vars, collapse = ", "), "\n",
      "Available variables include: ", paste(unique(df$variable), collapse = ", "),
      call. = FALSE
    )
  }
  
  # Keep optional name columns if they exist
  opt_name_cols <- intersect(names(df2), c("province_name", "district_name"))
  
  df_w <- df2 %>%
    dplyr::select(dplyr::all_of(c("district_id","year","month","variable","value", opt_name_cols))) %>%
    tidyr::pivot_wider(names_from = variable, values_from = value)
  
  temp_col <- intersect(names(df_w), temp_vars)[1]
  rh_col   <- intersect(names(df_w), rh_vars)[1]
  
  if (is.na(temp_col) || is.na(rh_col)) {
    stop(
      "After pivot, temp or RH column missing.\n",
      "Columns present: ", paste(names(df_w), collapse = ", "),
      call. = FALSE
    )
  }
  
  # Create standardized columns
  out <- df_w %>%
    dplyr::mutate(
      temp_raw = as.numeric(.data[[temp_col]]),
      rh = as.numeric(.data[[rh_col]])
    )
  
  # Kelvin heuristic
  if (stats::median(out$temp_raw, na.rm = TRUE) > 100) {
    message("ℹ Climate temperature appears to be Kelvin; converting to °C (T - 273.15).")
    out$temp_raw <- out$temp_raw - 273.15
  }
  
  # Deterministic, size-safe optional columns
  province_vec <- if ("province_name" %in% names(out)) as.character(out$province_name) else rep(NA_character_, nrow(out))
  district_vec <- if ("district_name" %in% names(out)) as.character(out$district_name) else rep(NA_character_, nrow(out))
  
  out <- out %>%
    dplyr::transmute(
      district_id = .data$district_id,
      year = .data$year,
      month = .data$month,
      temp_c = .data$temp_raw,
      rh = .data$rh,
      province_name = province_vec,
      district_name = district_vec
    )
  
  if (any(!is.finite(out$temp_c))) stop("df_climate has NA/non-finite temp_c.", call.=FALSE)
  if (any(!is.finite(out$rh))) stop("df_climate has NA/non-finite rh.", call.=FALSE)
  if (any(out$rh < 0 | out$rh > 100, na.rm = TRUE)) stop("df_climate rh outside 0..100.", call.=FALSE)
  
  assert_unique_key(out, c("district_id","year","month"), "df_climate")
  out
}

## ----------------------------------------------------------
## 5) Species mapping: prefer sentinel_species.csv else fallback IDs
## ----------------------------------------------------------
load_or_build_species_map <- function(district_ids) {
  if (length(district_ids) == 0) stop("No district_ids provided.", call.=FALSE)
  
  fp <- file.path(DIR_PROCESSED, "sentinel_species.csv")
  
  if (file.exists(fp)) {
    sp <- readr::read_csv(fp, show_col_types = FALSE)
    assert_has_cols(sp, c("district_id","species"), "sentinel_species.csv")
    sp$species <- tolower(as.character(sp$species))
    if (any(!sp$species %in% c("aegypti","albopictus"))) {
      stop("sentinel_species.csv species must be one of: aegypti, albopictus", call.=FALSE)
    }
    if (any(duplicated(sp$district_id))) stop("sentinel_species.csv has duplicate district_id", call.=FALSE)
    
    missing_ids <- setdiff(unique(district_ids), unique(sp$district_id))
    if (length(missing_ids) > 0) {
      stop(
        "sentinel_species.csv does not cover all districts. Missing:\n- ",
        paste(missing_ids, collapse = "\n- "),
        call. = FALSE
      )
    }
    return(sp)
  }
  
  # Fallback deterministic rule using your verified IDs:
  ids_aeg <- c("TUR.10.4_1", "TUR.81.6_1")  # Artvin + Zonguldak
  sp <- tibble::tibble(
    district_id = unique(district_ids),
    species = ifelse(unique(district_ids) %in% ids_aeg, "aegypti", "albopictus")
  )
  readr::write_csv(sp, fp)
  message("ℹ sentinel_species.csv not found; created fallback mapping at: ", fp)
  sp
}

## ----------------------------------------------------------
## 6) Core monthly CTMC computation
## ----------------------------------------------------------
## ----------------------------------------------------------
## 6) Core monthly CTMC computation
##    CORRECTED: beta_hv parameter added
## ----------------------------------------------------------
compute_ctmc_month_spark <- function(temp_c, rh, lambda_import_month,
                               tab_traits, beta_vh, beta_hv, m,
                               gamma_day, tau, i0 = 1L,
                               include_import_in_q = FALSE,
                               use_stochastic_EIP = TRUE,
                               eip_mc_n = 2000L,
                               eip_seed = NULL,
                               major_rule = list(type = "establishment")) {
  
  major_rule <- .major_rule_validate(major_rule)
  
  # Upfront guard: gamma_day must be strictly positive (avoids log(0) in extinction_prob_bd).
  if (!is.finite(gamma_day) || gamma_day <= 0) {
    stop("gamma_day must be a finite positive number (1 / infectious_period_days).", call. = FALSE)
  }
  
  # Double-counting guard: if include_import_in_q = TRUE, importation is already
  # folded into the BD birth rates, so multiplying by q_import_month afterwards
  # would count importation twice.  Emit a one-time warning to make the choice explicit.
  if (isTRUE(include_import_in_q) && lambda_import_month > 0) {
    warning(
      "include_import_in_q = TRUE: importation is included in the BD birth rates. ",
      "Do NOT also multiply p_local_major by q_import_month in the caller, ",
      "or importation will be double-counted.",
      call. = FALSE
    )
  }
  n <- 1:(tau - 1)
  mu_n <- gamma_day * n
  
  # ---------------------------------------------------------------
  # Jensen-correct EIP Monte Carlo (Appendix 1, Section 3)
  #
  # When use_stochastic_EIP = TRUE:
  #   1) Draw n_mc individual EIP values from LogNormal(T).
  #   2) For each draw j, compute λ_j = f(EIP_j) and then
  #      P_est_j via extinction_prob_bd.
  #   3) Return E[P_est] = mean(P_est_j).
  #
  # Jensen direction depends on R_eff regime:
  #   R_eff >> 1: P_est concavity dominates → E[P_est] ≤ P_est(E[EIP])
  #              (conservative / lower)
  #   R_eff ≈  1: max(0,...) threshold convexity dominates →
  #              E[P_est] ≥ P_est(E[EIP])
  # In both cases |bias| < 5% for σ_log = 0.45.
  #
  # Previous code averaged survival E[exp(-μv·EIP)] first, then
  # computed a single P_est.  This version averages at the P_est
  # level, which is conceptually correct for individual EIP
  # heterogeneity in the stochastic spark phase.
  # ---------------------------------------------------------------
  
  # Build base call_args (shared between draw-level and scalar paths)
  call_args_base <- list(
    n_infected = 1L,
    T = temp_c,
    RH = rh,
    tab_traits = tab_traits,
    beta_vh = beta_vh,
    beta_hv = beta_hv,
    m = m
  )
  
  if (isTRUE(use_stochastic_EIP)) {
    # --- Stochastic EIP: draw-level P_est averaging ---
    call_args <- c(call_args_base, list(
      use_stochastic_EIP = TRUE,
      n_mc = eip_mc_n,
      seed = eip_seed,
      return_draws = TRUE      # ← get individual λ per EIP draw
    ))
    
    lam_draws <- do.call(make_lambda_local, call_args)
    if (any(!is.finite(lam_draws))) {
      lam_draws[!is.finite(lam_draws)] <- 0
    }
    
    # Compute P_est for EACH EIP draw, then average (Jensen-correct)
    .compute_p_est_one <- function(li) {
      if (li <= 0) return(0)
      lambda_local_n <- n * li
      lambda_birth_n <- if (isTRUE(include_import_in_q)) {
        lambda_local_n + lambda_import_month
      } else {
        lambda_local_n
      }
      qi <- extinction_prob_bd(lambda_n = lambda_birth_n, mu_n = mu_n,
                               tau = tau, i0 = i0)
      1 - qi
    }
    
    p_est_draws <- vapply(lam_draws, .compute_p_est_one, numeric(1))
    p_est   <- mean(p_est_draws)       # E[P_est(EIP_i)]
    q       <- 1 - p_est
    lam_i1  <- mean(lam_draws)         # diagnostic: E[λ(EIP_i)]
    
  } else {
    # --- Deterministic EIP: single scalar λ ---
    call_args <- c(call_args_base, list(
      use_stochastic_EIP = FALSE
    ))
    
    lam_i1 <- do.call(make_lambda_local, call_args)
    if (!is.finite(lam_i1) || is.na(lam_i1)) {
      stop("make_lambda_local returned NA/non-finite for i=1.", call. = FALSE)
    }
    if (lam_i1 < 0) {
      stop("make_lambda_local returned negative value for i=1.", call. = FALSE)
    }
    
    lambda_local_n <- n * lam_i1
    lambda_birth_n <- if (isTRUE(include_import_in_q)) {
      lambda_local_n + lambda_import_month
    } else {
      lambda_local_n
    }
    
    q     <- extinction_prob_bd(lambda_n = lambda_birth_n, mu_n = mu_n,
                                tau = tau, i0 = i0)
    p_est <- 1 - q
  }
  
  # --- Local major probability (import OFF) ---
  # For stochastic EIP + non-establishment major rules, use mean λ
  # (Gillespie paths already incorporate stochasticity internally).
  p_local_major <- switch(
    major_rule$type,
    "establishment" = p_est,
    "hit_K" = {
      K <- as.integer(major_rule$K)
      nK <- 1:(K - 1L)
      qK <- extinction_prob_bd(lambda_n = nK * lam_i1, mu_n = gamma_day * nK, tau = K, i0 = 1L)
      1 - qK
    },
    "sustain_days" = {
      .sim_bd_prob(lambda_local_i1 = lam_i1, gamma_day = gamma_day,
                   D_target = as.numeric(major_rule$D_days), X_target = NULL,
                   n_path = as.integer(major_rule$n_path),
                   seed = if (!is.null(eip_seed)) as.integer(eip_seed + 101L) else NULL)
    },
    "secondary_ge_X" = {
      .sim_bd_prob(lambda_local_i1 = lam_i1, gamma_day = gamma_day,
                   D_target = NULL, X_target = as.integer(major_rule$X),
                   n_path = as.integer(major_rule$n_path),
                   seed = if (!is.null(eip_seed)) as.integer(eip_seed + 202L) else NULL)
    },
    stop("unknown major_rule$type")
  )
  
  list(
    p_local_major   = p_local_major,
    p_establishment = p_est,
    q_extinction    = q,
    lambda_local_i1 = lam_i1
  )
}
## ----------------------------------------------------------
## 7) Main runner
##    CORRECTED: beta_hv parameter added
## ----------------------------------------------------------
run_ctmc_spark <- function(
    year_min = 2025,
    year_max = 2075,
    tau = 30L,
    infectious_period_days = 5,
    beta_vh = 0.3,
    beta_hv = 0.33,  # ← CORRECTED: beta_hv parameter added with default value
    m = 1.0,         # Ross-Macdonald standard (1 female mosquito per person)
    include_import_in_q = FALSE,
    # ---- NEW: stochastic EIP toggle (single control point)
    use_stochastic_EIP = TRUE,   # tezde log-normal EIP kullanılmıştır
    eip_mc_n = 2000L,
    eip_seed = NULL,
    major_rule = list(type = "establishment"),   # NEW
    out_dir = file.path(DIR_OUTPUT_SSP, "model_results")
) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fp_import  <- file.path(DIR_PROCESSED_SSP, "importation_pressure_monthly_2025_2075.rds")
  fp_climate <- file.path(DIR_PROCESSED_SSP, "climate_sentinel_monthly_long_2015_2100.rds")
  
  fp_traits_aeg <- file.path(DIR_PROCESSED, "trait_params_aegypti.csv")
  fp_traits_alb <- file.path(DIR_PROCESSED, "trait_params_albopictus.csv")
  
  if (!file.exists(fp_import))  stop("Missing: ", fp_import,  call.=FALSE)
  if (!file.exists(fp_climate)) stop("Missing: ", fp_climate, call.=FALSE)
  if (!file.exists(fp_traits_aeg)) stop("Missing: ", fp_traits_aeg, call.=FALSE)
  if (!file.exists(fp_traits_alb)) stop("Missing: ", fp_traits_alb, call.=FALSE)
  
  message("ℹ Inputs:")
  message("- importation: ", fp_import)
  message("- climate:     ", fp_climate)
  message("- traits aeg:  ", fp_traits_aeg)
  message("- traits alb:  ", fp_traits_alb)
  
  df_import  <- load_importation_monthly(fp_import)
  df_climate <- load_climate_monthly_long(fp_climate)
  
  # ---- Sanity check (INSIDE function!)
  message("ℹ Sanity check:")
  message("  - Climate districts: ", dplyr::n_distinct(df_climate$district_id))
  message("  - Import districts:  ", dplyr::n_distinct(df_import$district_id))
  message("  - Overlap districts: ",
          length(intersect(unique(df_climate$district_id), unique(df_import$district_id))))
  
  # ---- Filter time window
  df_import  <- df_import  %>% dplyr::filter(.data$year >= year_min, .data$year <= year_max)
  df_climate <- df_climate %>% dplyr::filter(.data$year >= year_min, .data$year <= year_max)
  
  # ---- Join
  key <- c("district_id","year","month")
  df <- df_climate %>%
    dplyr::left_join(df_import, by = key)
  
  if (any(is.na(df$lambda_import))) {
    message("ℹ Missing lambda_import for some rows; setting missing to 0.")
    df$lambda_import[is.na(df$lambda_import)] <- 0
  }
  
  if (any(is.na(df$q_import_month))) {
    message("ℹ Missing q_import_month for some rows; setting missing to 0.")
    df$q_import_month[is.na(df$q_import_month)] <- 0
  }
  
  # ---- Species mapping
  sp_map <- load_or_build_species_map(unique(df$district_id))
  df <- df %>% dplyr::left_join(sp_map, by = "district_id")
  if (any(is.na(df$species))) stop("species is NA after join; check sentinel_species mapping.", call.=FALSE)
  
  # ---- Trait tables
  tab_traits_aeg <- read_trait_params(fp_traits_aeg)
  tab_traits_alb <- read_trait_params(fp_traits_alb)
  
  gamma_day <- 1 / infectious_period_days
  
  # ---- Monthly CTMC
  res <- df %>%
    dplyr::mutate(
      tmp = purrr::pmap(
        list(district_id, year, month, temp_c, rh, lambda_import, species),
        function(district_id, year, month, temp_c, rh, lambda_import, species) {
          tab_traits <- if (species == "aegypti") tab_traits_aeg else tab_traits_alb
          
          seed_here <- NULL
          if (isTRUE(use_stochastic_EIP) && !is.null(eip_seed)) {
            # district_id + year + month üzerinden deterministik seed üret
            key_str <- paste(district_id, year, month, sep = "|")
            seed_here <- as.integer((abs(sum(utf8ToInt(key_str))) + as.integer(eip_seed)) %% .Machine$integer.max)
            if (seed_here == 0L) seed_here <- 1L
          }
          
          compute_ctmc_month_spark(
            temp_c = temp_c, rh = rh, lambda_import_month = lambda_import,
            tab_traits = tab_traits, 
            beta_vh = beta_vh,
            beta_hv = beta_hv,  # ← CORRECTED: beta_hv parameter passed
            m = m,
            gamma_day = gamma_day, tau = tau, i0 = 1L,
            include_import_in_q = include_import_in_q,
            use_stochastic_EIP = use_stochastic_EIP,
            eip_mc_n = eip_mc_n, eip_seed = seed_here,
            major_rule = major_rule
          )
        }
      ),
      p_local_major   = purrr::map_dbl(tmp, "p_local_major"),   # NEW
      p_establishment = purrr::map_dbl(tmp, "p_establishment"),
      q_extinction    = purrr::map_dbl(tmp, "q_extinction"),
      lambda_local_i1 = purrr::map_dbl(tmp, "lambda_local_i1")
    ) %>%
    dplyr::select(-tmp) %>%
    dplyr::mutate(
      # Aylık majör olasılık = ithal >=1 olasılığı * yerel major olasılığı
      p_month_major = pmin(pmax(.data$q_import_month, 0), 1) * pmin(pmax(.data$p_local_major, 0), 1)
    )
  
  # ---- Yearly p(>=1 major outbreak) using Bernoulli months:
  # p_year = 1 - prod_m (1 - p_month_major)
  res_yearly <- res %>%
    dplyr::group_by(.data$district_id, .data$year) %>%
    dplyr::summarise(
      species = dplyr::first(.data$species),
      
      # Trace/QA
      Lambda_import_year  = sum(.data$lambda_import, na.rm = TRUE),
      mean_p_est_year     = mean(.data$p_establishment, na.rm = TRUE),
      mean_p_local_major  = mean(.data$p_local_major, na.rm = TRUE),
      
      # Bernoulli aggregation
      p_ge1_major_year    = 1 - prod(1 - pmin(pmax(.data$p_month_major, 0), 1), na.rm = TRUE),
      
      # Diagnostics
      mean_q_import_year  = mean(.data$q_import_month, na.rm = TRUE),
      mean_p_month_major  = mean(.data$p_month_major, na.rm = TRUE),
      
      mean_temp_c_year    = mean(.data$temp_c, na.rm = TRUE),
      mean_rh_year        = mean(.data$rh, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- Whole horizon per sentinel
  res_horizon <- res %>%
    dplyr::group_by(.data$district_id) %>%
    dplyr::summarise(
      species = dplyr::first(.data$species),
      year_min = min(.data$year),
      year_max = max(.data$year),
      
      Lambda_import     = sum(.data$lambda_import, na.rm = TRUE),
      p_ge1_major       = 1 - prod(1 - pmin(pmax(.data$p_month_major, 0), 1), na.rm = TRUE),
      
      mean_p_est        = mean(.data$p_establishment, na.rm = TRUE),
      mean_p_local_major= mean(.data$p_local_major,   na.rm = TRUE),
      mean_q_import     = mean(.data$q_import_month,  na.rm = TRUE),
      mean_p_month_major= mean(.data$p_month_major,   na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- Write outputs
  fp_monthly_csv <- file.path(out_dir, sprintf("ctmc_spark_monthly_%d_%d.csv", year_min, year_max))
  fp_monthly_rds <- file.path(out_dir, sprintf("ctmc_spark_monthly_%d_%d.rds", year_min, year_max))
  fp_yearly_csv  <- file.path(out_dir, sprintf("ctmc_spark_yearly_%d_%d.csv",  year_min, year_max))
  fp_yearly_rds  <- file.path(out_dir, sprintf("ctmc_spark_yearly_%d_%d.rds",  year_min, year_max))
  fp_horizon_csv <- file.path(out_dir, sprintf("ctmc_spark_horizon_%d_%d.csv", year_min, year_max))
  fp_horizon_rds <- file.path(out_dir, sprintf("ctmc_spark_horizon_%d_%d.rds", year_min, year_max))
  
  readr::write_csv(res,         fp_monthly_csv)
  saveRDS(res,                  fp_monthly_rds)
  readr::write_csv(res_yearly,  fp_yearly_csv)
  saveRDS(res_yearly,           fp_yearly_rds)
  readr::write_csv(res_horizon, fp_horizon_csv)
  saveRDS(res_horizon,          fp_horizon_rds)
  
  message("✅ Outputs written:")
  message("- monthly: ",  fp_monthly_csv)
  message("- yearly:  ",  fp_yearly_csv)
  message("- horizon: ",  fp_horizon_csv)
  
  invisible(list(monthly = res, yearly = res_yearly, horizon = res_horizon))
}

## ----------------------------------------------------------
## 8) Entrypoint
##    CORRECTED: beta_hv and m parameters updated
##    m = 1.0 (Ross-Macdonald standard); sensitivity at {0.5, 0.8, 1.2, 2.0}
## ----------------------------------------------------------
if (identical(environment(), globalenv())) {
  run_ctmc_spark(
    year_min = 2025,
    year_max = 2075,
    tau = 30L,
    infectious_period_days = 5,
    beta_vh = 0.3,
    beta_hv = 0.33,  # ← CORRECTED: beta_hv parameter added
    m = 1.0,         # Ross-Macdonald standard (Section 3.5, Table 1 — UPDATED)
    include_import_in_q = FALSE,
    use_stochastic_EIP = TRUE,
    eip_mc_n = 2000L,
    eip_seed = 123,
    # Örnek kullanım (analitik eşik): I ≥ 20'ye ulaşırsa major
    # major_rule = list(type = "hit_K", K = 20),
    out_dir = file.path(DIR_OUTPUT_SSP, "model_results")
  )
}

message("✓ ctmc_spark.R loaded (SSP-aware paths + gamma_day assert + double-counting guard)")