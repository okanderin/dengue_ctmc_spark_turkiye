# ==========================================================
# R/03_models/ctmc_spark_monte_carlo.R
# Spark-phase CTMC (birth-death) + species-aware traits + MC replication
# - Climate: long format (variable/value) -> temp_c, rh
# - Importation: monthly expected imports (lambda_import)
# - Species rule: sentinel_species.csv (preferred) / fallback deterministic
# - Stochastic EIP (log-normal) + n_rep Monte Carlo outer loop
# - Outputs: monthly + yearly + horizon risk with mean/p2.5/p97.5
# - Filenames suffix: _rep{n_rep}
# - SSP-AWARE: reads from DIR_PROCESSED_SSP, writes to DIR_OUTPUT_SSP
#
# CORRECTED VERSION: beta_hv parameter added to lambda calculation
# ==========================================================

## ----------------------------------------------------------
## 0) Initialize project environment
## ----------------------------------------------------------
source("R/01_setup/init.R")
source("R/03_models/parameter_functions.R")  # must define: read_trait_params(), make_lambda_local()

# --- Ensure importation pressure path exists (SSP-aware) ---
if (!exists("fp_import")) {
  fp_import_guess <- file.path(DIR_PROCESSED_SSP, "importation_pressure_monthly_2025_2075.rds")
  if (!file.exists(fp_import_guess)) {
    # Fallback: try legacy filename
    fp_import_guess2 <- file.path(DIR_PROCESSED_SSP, "importation_pressure_monthly_2024_2100.rds")
    if (file.exists(fp_import_guess2)) {
      warning("Using legacy importation file: ", fp_import_guess2,
              "\nConsider re-running build_importation_pressure_monthly.R (v3).",
              call. = FALSE)
      fp_import_guess <- fp_import_guess2
    } else {
      stop(
        "fp_import is not defined and no importation file found:\n",
        "  Tried: ", fp_import_guess, "\n",
        "  Tried: ", fp_import_guess2, "\n",
        "Fix: run build_importation_pressure_monthly.R for this SSP first.",
        call. = FALSE
      )
    }
  }
  fp_import <- fp_import_guess
}

## Deterministic row-level seed helper (no extra packages)
.hash_seed <- function(did, year, month, rep_id, base = 104729L) {
  s <- sum(utf8ToInt(as.character(did))) +
    as.integer(year) * 1009L + as.integer(month) * 29L +
    as.integer(rep_id) * 7919L
  s <- as.integer(abs(s %% .Machine$integer.max))
  if (s == 0L) s <- base
  s
}

## ----------------------------------------------------------
## major_rule validator (mirrors ctmc_spark.R exactly)
## type = "establishment" | "hit_K" | "sustain_days" | "secondary_ge_X"
## ----------------------------------------------------------
.major_rule_validate <- function(major_rule) {
  if (is.null(major_rule) || is.null(major_rule$type)) {
    return(list(type = "establishment"))
  }
  t <- match.arg(major_rule$type,
                 c("establishment", "hit_K", "sustain_days", "secondary_ge_X"))
  out <- list(type = t)
  if (t == "hit_K") {
    K <- as.integer(major_rule$K %||% 20L)
    if (is.na(K) || K < 2L) stop("major_rule$K must be >= 2 for hit_K.", call. = FALSE)
    out$K <- K
  } else if (t == "sustain_days") {
    D      <- as.numeric(major_rule$D_days  %||% 14)
    n_path <- as.integer(major_rule$n_path  %||% 2000L)
    if (!is.finite(D) || D <= 0) stop("major_rule$D_days must be > 0.", call. = FALSE)
    if (n_path < 100L) warning("n_path is small; increase to >1000 for stability.", call. = FALSE)
    out$D_days <- D; out$n_path <- n_path; out$seed <- major_rule$seed %||% NULL
  } else if (t == "secondary_ge_X") {
    X      <- as.integer(major_rule$X      %||% 20L)
    n_path <- as.integer(major_rule$n_path %||% 2000L)
    if (is.na(X) || X < 1L) stop("major_rule$X must be >= 1.", call. = FALSE)
    out$X <- X; out$n_path <- n_path; out$seed <- major_rule$seed %||% NULL
  }
  out
}

## Null-coalescing helper (used by .major_rule_validate)
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ----------------------------------------------------------
## Gillespie BD single-path simulator (I(0)=1, import OFF)
## Used for sustain_days / secondary_ge_X major_rule types.
## ----------------------------------------------------------
.sim_bd_once <- function(lambda_local_i1, gamma_day,
                         D_target = NULL, X_target = NULL,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n         <- 1L
  t         <- 0.0
  sec_cases <- 0L
  while (n > 0L) {
    lam <- n * lambda_local_i1
    mu  <- n * gamma_day
    tot <- lam + mu
    if (tot <= 0) break
    dt <- rexp(1L, rate = tot)
    t  <- t + dt
    if (!is.null(D_target) && t >= D_target)
      return(list(hit = TRUE, t = t, sec = sec_cases))
    if (runif(1) < lam / tot) {
      n <- n + 1L
      sec_cases <- sec_cases + 1L
      if (!is.null(X_target) && sec_cases >= X_target)
        return(list(hit = TRUE, t = t, sec = sec_cases))
    } else {
      n <- n - 1L
    }
  }
  list(hit = FALSE, t = t, sec = sec_cases)
}

.sim_bd_prob <- function(lambda_local_i1, gamma_day,
                         D_target = NULL, X_target = NULL,
                         n_path = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  hits <- 0L
  for (r in seq_len(n_path)) {
    s_here <- if (is.null(seed)) NULL else
      as.integer((seed + r) %% .Machine$integer.max)
    res <- .sim_bd_once(lambda_local_i1, gamma_day,
                        D_target, X_target, seed = s_here)
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
## 2)Birth–death extinction probability (absorption at 0 vs tau)
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
##    (Harmonized with ctmc_spark.R)
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
    df$q_import_month <- 1 - exp(
      -as.numeric(df$lambda_import_per_day) *
        pmin(as.numeric(df$d_days), as.numeric(df$days_in_month))
    )
  } else if ("lambda_import" %in% names(df) && "days_in_month" %in% names(df)) {
    # fallback: interpret lambda_import as expected/month and approximate daily by /days (3-day window)
    df$q_import_month <- 1 - exp(
      -(as.numeric(df$lambda_import) / as.numeric(df$days_in_month)) *
        pmin(3, as.numeric(df$days_in_month))
    )
  } else {
    # as a last resort, set to 0
    df$q_import_month <- 0
  }
  assert_numeric_nonneg(df$q_import_month, "df_import$q_import_month")
  df$q_import_month <- pmin(pmax(df$q_import_month, 0), 1)
  
  # days_in_month: if missing or NA, compute
  if (!("days_in_month" %in% names(df))) {
    df$days_in_month <- as.integer(lubridate::days_in_month(
      lubridate::make_date(df$year, df$month, 1)
    ))
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
    keep <- c(
      key, "lambda_import", "q_import_month", "days_in_month",
      intersect(names(df), c("province_name", "district_name"))
    )
    df <- df %>% dplyr::select(dplyr::all_of(keep))
  }
  
  assert_unique_key(df, key, "df_import")
  df
}

## ----------------------------------------------------------
## 4) Load climate (monthly long) -> temp_c + rh
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
  
  out <- df_w %>%
    dplyr::mutate(
      temp_raw = as.numeric(.data[[temp_col]]),
      rh = as.numeric(.data[[rh_col]])
    )
  
  if (stats::median(out$temp_raw, na.rm = TRUE) > 100) {
    message("ℹ Climate temperature appears to be Kelvin; converting to °C (T - 273.15).")
    out$temp_raw <- out$temp_raw - 273.15
  }
  
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
  if (any(!is.finite(out$rh)))    stop("df_climate has NA/non-finite rh.", call.=FALSE)
  if (any(out$rh < 0 | out$rh > 100, na.rm = TRUE)) stop("df_climate rh outside 0..100.", call.=FALSE)
  
  assert_unique_key(out, c("district_id","year","month"), "df_climate")
  out
}

## ----------------------------------------------------------
## 5) Species mapping
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
  
  ids_aeg <- c("TUR.10.4_1", "TUR.81.6_1")  # Artvin + Zonguldak
  ids_alb <- setdiff(unique(district_ids), ids_aeg)
  
  sp <- data.frame(
    district_id = c(ids_aeg, ids_alb),
    species = c(rep("aegypti", length(ids_aeg)), rep("albopictus", length(ids_alb))),
    stringsAsFactors = FALSE
  )
  
  message("⚠ sentinel_species.csv not found; using fallback deterministic rule:")
  message("  - aegypti:   ", paste(ids_aeg, collapse=", "))
  message("  - albopictus: ", paste(ids_alb, collapse=", "))
  sp
}


## ----------------------------------------------------------
## 6) Core monthly CTMC computation (ONE row; stochastic EIP enabled)
##    CORRECTED: beta_hv parameter added
##    UPDATED:   major_rule support + p_local_major (mirrors ctmc_spark.R)
## ----------------------------------------------------------
compute_ctmc_month_mc <- function(district_id, year, month,
                               temp_c, rh, lambda_import_month,
                               tab_traits, beta_vh, beta_hv, m,
                               gamma_day, tau, i0 = 1L,
                               include_import_in_q = FALSE,
                               use_stochastic_EIP = TRUE,
                               n_mc_eip = 2000L,
                               rep_id = 1L,
                               base_seed = 123L,
                               major_rule = list(type = "establishment")) {
  
  major_rule <- .major_rule_validate(major_rule)
  
  # Upfront guard: gamma_day must be strictly positive.
  if (!is.finite(gamma_day) || gamma_day <= 0) {
    stop("gamma_day must be a finite positive number (1 / infectious_period_days).", call. = FALSE)
  }
  
  # Double-counting guard (same logic as ctmc_spark.R).
  if (isTRUE(include_import_in_q) && lambda_import_month > 0) {
    warning(
      "include_import_in_q = TRUE: importation is included in the BD birth rates. ",
      "Do NOT also multiply p_local_major by q_import_month in the caller.",
      call. = FALSE
    )
  }
  
  n <- 1:(tau - 1)
  
  row_seed <- .hash_seed(district_id, year, month, rep_id)
  
  lam_i1 <- make_lambda_local(
    n_infected = 1L,
    T = temp_c,
    RH = rh,
    tab_traits = tab_traits,
    beta_vh = beta_vh,
    beta_hv = beta_hv,
    m = m,
    use_stochastic_EIP = use_stochastic_EIP,
    n_mc = n_mc_eip,
    seed = as.integer(base_seed + row_seed)
  )
  
  if (!is.finite(lam_i1) || is.na(lam_i1)) stop("make_lambda_local returned NA/non-finite for i=1.", call.=FALSE)
  if (lam_i1 < 0) stop("make_lambda_local returned negative value for i=1.", call.=FALSE)
  
  lambda_local_n <- n * lam_i1
  lambda_birth_n <- if (isTRUE(include_import_in_q)) lambda_local_n + lambda_import_month else lambda_local_n
  mu_n <- gamma_day * n
  
  q     <- extinction_prob_bd(lambda_n = lambda_birth_n, mu_n = mu_n, tau = tau, i0 = i0)
  p_est <- 1 - q
  
  # --- p_local_major: respects major_rule (mirrors ctmc_spark.R) ---
  p_local_major <- switch(
    major_rule$type,
    "establishment" = p_est,
    "hit_K" = {
      K  <- as.integer(major_rule$K)
      nK <- 1:(K - 1L)
      qK <- extinction_prob_bd(lambda_n = nK * lam_i1, mu_n = gamma_day * nK,
                               tau = K, i0 = 1L)
      1 - qK
    },
    "sustain_days" = {
      .sim_bd_prob(lambda_local_i1 = lam_i1, gamma_day = gamma_day,
                   D_target = as.numeric(major_rule$D_days), X_target = NULL,
                   n_path   = as.integer(major_rule$n_path),
                   seed = if (!is.null(row_seed)) as.integer(row_seed + 101L) else NULL)
    },
    "secondary_ge_X" = {
      .sim_bd_prob(lambda_local_i1 = lam_i1, gamma_day = gamma_day,
                   D_target = NULL, X_target = as.integer(major_rule$X),
                   n_path   = as.integer(major_rule$n_path),
                   seed = if (!is.null(row_seed)) as.integer(row_seed + 202L) else NULL)
    },
    stop("unknown major_rule$type", call. = FALSE)
  )
  
  list(
    p_local_major   = p_local_major,   # used for p_month_major in caller
    p_establishment = p_est,           # kept for diagnostics / backward compat
    q_extinction    = q,
    lambda_local_i1 = lam_i1
  )
}

## ----------------------------------------------------------
## 7) Main runner with Monte Carlo replication
##    CORRECTED: beta_hv parameter added
##    UPDATED:   m = 1.0 (Ross-Macdonald standard); sensitivity at {0.5, 0.8, 1.2, 2.0}
## ----------------------------------------------------------
run_ctmc_spark <- function(
    year_min = 2025,
    year_max = 2075,
    tau = 30L,
    infectious_period_days = 5,
    beta_vh = 0.3,
    beta_hv = 0.33,
    m = 1.0,                              # Ross-Macdonald standard (1 female mosquito per person)
    include_import_in_q = FALSE,
    use_stochastic_EIP = TRUE,
    n_mc_eip = 2000L,                     # inner EIP MC
    n_rep = 1000L,                        # outer MC
    seed = 123L,                          # base seed (outer)
    major_rule = list(type = "establishment"),
    out_dir = file.path(DIR_OUTPUT_SSP, "simulation")
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # SSP-dependent inputs
  # NOTE: paths.R defines fp_import globally; use it if available,

  #       otherwise fall back to the canonical v3 filename (_2025_2075).
  #       The old _2024_2100 filename is from pre-GBD pipeline and should

  #       NOT be used for new runs.
  fp_import  <- file.path(DIR_PROCESSED_SSP, "importation_pressure_monthly_2025_2075.rds")
  fp_climate <- file.path(DIR_PROCESSED_SSP, "climate_sentinel_monthly_long_2015_2100.rds")
  
  # SSP-independent inputs (shared across scenarios)
  fp_traits_aeg <- file.path(DIR_PROCESSED, "trait_params_aegypti.csv")
  fp_traits_alb <- file.path(DIR_PROCESSED, "trait_params_albopictus.csv")
  
  if (!file.exists(fp_import)) stop("Missing: ", fp_import, call.=FALSE)
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
  
  message("ℹ Sanity check:")
  message("  - Climate districts: ", dplyr::n_distinct(df_climate$district_id))
  message("  - Import districts:  ", dplyr::n_distinct(df_import$district_id))
  message("  - Overlap districts: ",
          length(intersect(unique(df_climate$district_id), unique(df_import$district_id))))
  
  df_import  <- df_import  %>% dplyr::filter(.data$year >= year_min, .data$year <= year_max)
  df_climate <- df_climate %>% dplyr::filter(.data$year >= year_min, .data$year <= year_max)
  
  key <- c("district_id","year","month")
  df <- df_climate %>% dplyr::left_join(df_import, by = key)
  
  if (any(is.na(df$lambda_import))) {
    message("ℹ Missing lambda_import for some rows; setting missing to 0.")
    df$lambda_import[is.na(df$lambda_import)] <- 0
  }
  
  if (any(is.na(df$q_import_month))) {
    message("ℹ Missing q_import_month for some rows; setting missing to 0.")
    df$q_import_month[is.na(df$q_import_month)] <- 0
  }
  
  sp_map <- load_or_build_species_map(unique(df$district_id))
  df <- df %>% dplyr::left_join(sp_map, by = "district_id")
  if (any(is.na(df$species))) stop("species is NA after join; check sentinel_species mapping.", call.=FALSE)
  
  tab_traits_aeg <- read_trait_params(fp_traits_aeg)
  tab_traits_alb <- read_trait_params(fp_traits_alb)
  
  gamma_day <- 1 / infectious_period_days
  
  set.seed(seed)
  rep_ids <- seq_len(as.integer(n_rep))
  
  compute_one_rep <- function(rep_id) {
    res_month <- df %>%
      dplyr::mutate(
        tmp = purrr::pmap(
          list(district_id, year, month, temp_c, rh, lambda_import, species),
          function(district_id, year, month, temp_c, rh, lambda_import, species) {
            tab_traits <- if (species == "aegypti") tab_traits_aeg else tab_traits_alb
            compute_ctmc_month_mc(
              district_id = district_id,
              year = year, month = month,
              temp_c = temp_c,
              rh = rh,
              lambda_import_month = lambda_import,
              tab_traits = tab_traits,
              beta_vh = beta_vh,
              beta_hv = beta_hv,
              m = m,
              gamma_day = gamma_day,
              tau = tau,
              i0 = 1L,
              include_import_in_q = include_import_in_q,
              use_stochastic_EIP = use_stochastic_EIP,
              n_mc_eip = n_mc_eip,
              rep_id = rep_id,
              base_seed = seed,
              major_rule = major_rule   # passed through from run_ctmc_spark()
            )
          }
        ),
        p_local_major   = purrr::map_dbl(tmp, "p_local_major"),   # respects major_rule
        p_establishment = purrr::map_dbl(tmp, "p_establishment"),  # kept for diagnostics
        q_extinction    = purrr::map_dbl(tmp, "q_extinction"),
        lambda_local_i1 = purrr::map_dbl(tmp, "lambda_local_i1"),
        # p_month_major now uses p_local_major (not p_establishment directly),
        # consistent with ctmc_spark.R and respecting alternative major_rule types.
        p_month_major = pmin(pmax(.data$q_import_month, 0), 1) *
          pmin(pmax(.data$p_local_major, 0), 1)
      ) %>%
      dplyr::select(-tmp) %>%
      dplyr::mutate(.rep = rep_id)
    
    res_year <- res_month %>%
      dplyr::group_by(.data$district_id, .data$year) %>%
      dplyr::summarise(
        species = dplyr::first(.data$species),
        
        # trace
        Lambda_import_year = sum(.data$lambda_import, na.rm = TRUE),
        
        # Bernoulli aggregation
        p_ge1_major_year = 1 - prod(1 - pmin(pmax(.data$p_month_major, 0), 1), na.rm = TRUE),
        
        # diagnostics
        mean_q_import_year   = mean(.data$q_import_month, na.rm = TRUE),
        mean_p_local_major   = mean(.data$p_local_major,   na.rm = TRUE),  # NEW
        mean_p_est_year      = mean(.data$p_establishment, na.rm = TRUE),
        mean_p_month_major   = mean(.data$p_month_major, na.rm = TRUE),
        
        mean_temp_c_year = mean(.data$temp_c, na.rm = TRUE),
        mean_rh_year     = mean(.data$rh, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(.rep = rep_id)
    
    res_hor <- res_month %>%
      dplyr::group_by(.data$district_id) %>%
      dplyr::summarise(
        species  = dplyr::first(.data$species),
        year_min = min(.data$year),
        year_max = max(.data$year),
        
        Lambda_import = sum(.data$lambda_import, na.rm = TRUE),
        
        p_ge1_major = 1 - prod(1 - pmin(pmax(.data$p_month_major, 0), 1), na.rm = TRUE),
        
        mean_q_import      = mean(.data$q_import_month, na.rm = TRUE),
        mean_p_local_major = mean(.data$p_local_major,   na.rm = TRUE),  # NEW
        mean_p_est         = mean(.data$p_establishment, na.rm = TRUE),
        mean_p_month_major = mean(.data$p_month_major, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(.rep = rep_id)
    
    list(monthly = res_month, yearly = res_year, horizon = res_hor)
  }
  
  message("ℹ Running Monte Carlo reps: n_rep = ", n_rep, " (inner EIP MC: ", n_mc_eip, ")")
  reps <- lapply(rep_ids, compute_one_rep)
  
  monthly_all <- dplyr::bind_rows(lapply(reps, `[[`, "monthly"))
  yearly_all  <- dplyr::bind_rows(lapply(reps, `[[`, "yearly"))
  horizon_all <- dplyr::bind_rows(lapply(reps, `[[`, "horizon"))
  
  qfun <- function(x, p) stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE)
  
  monthly_sum <- monthly_all %>%
    dplyr::group_by(district_id, year, month) %>%
    dplyr::summarise(
      species = dplyr::first(species),
      temp_c  = dplyr::first(temp_c),
      rh      = dplyr::first(rh),
      lambda_import = dplyr::first(lambda_import),
      
      lambda_local_i1_mean  = mean(lambda_local_i1, na.rm = TRUE),
      lambda_local_i1_p2_5  = qfun(lambda_local_i1, 0.025),
      lambda_local_i1_p97_5 = qfun(lambda_local_i1, 0.975),
      
      p_establishment_mean  = mean(p_establishment, na.rm = TRUE),
      p_establishment_p2_5  = qfun(p_establishment, 0.025),
      p_establishment_p97_5 = qfun(p_establishment, 0.975),
      
      q_extinction_mean     = mean(q_extinction, na.rm = TRUE),
      q_extinction_p2_5     = qfun(q_extinction, 0.025),
      q_extinction_p97_5    = qfun(q_extinction, 0.975),
      
      q_import_month        = dplyr::first(q_import_month),
      p_month_major_mean    = mean(p_month_major, na.rm = TRUE),
      p_month_major_p2_5    = qfun(p_month_major, 0.025),
      p_month_major_p97_5   = qfun(p_month_major, 0.975),
      
      .groups = "drop"
    )
  
  qa_bad <- monthly_sum %>%
    dplyr::filter(q_import_month < 0 | q_import_month > 1 | !is.finite(q_import_month))
  if (nrow(qa_bad) > 0) stop("QA failed: q_import_month outside [0,1].", call. = FALSE)
  
  yearly_sum <- yearly_all %>%
    dplyr::group_by(district_id, year) %>%
    dplyr::summarise(
      species = dplyr::first(species),
      Lambda_import_year = dplyr::first(Lambda_import_year),
      
      p_ge1_major_year_mean  = mean(p_ge1_major_year, na.rm = TRUE),
      p_ge1_major_year_p2_5  = qfun(p_ge1_major_year, 0.025),
      p_ge1_major_year_p97_5 = qfun(p_ge1_major_year, 0.975),
      
      mean_q_import_year_mean = mean(mean_q_import_year, na.rm = TRUE),
      mean_q_import_year_p2_5 = qfun(mean_q_import_year, 0.025),
      mean_q_import_year_p97_5= qfun(mean_q_import_year, 0.975),
      
      mean_p_month_major_mean = mean(mean_p_month_major, na.rm = TRUE),
      mean_p_month_major_p2_5 = qfun(mean_p_month_major, 0.025),
      mean_p_month_major_p97_5= qfun(mean_p_month_major, 0.975),
      
      mean_p_est_year_mean    = mean(mean_p_est_year, na.rm = TRUE),
      mean_p_est_year_p2_5    = qfun(mean_p_est_year, 0.025),
      mean_p_est_year_p97_5   = qfun(mean_p_est_year, 0.975),
      
      mean_temp_c_year = dplyr::first(mean_temp_c_year),
      mean_rh_year     = dplyr::first(mean_rh_year),
      .groups = "drop"
    )
  
  horizon_sum <- horizon_all %>%
    dplyr::group_by(district_id) %>%
    dplyr::summarise(
      species   = dplyr::first(species),
      year_min  = dplyr::first(year_min),
      year_max  = dplyr::first(year_max),
      Lambda_import = dplyr::first(Lambda_import),
      
      p_ge1_major_mean  = mean(p_ge1_major, na.rm = TRUE),
      p_ge1_major_p2_5  = qfun(p_ge1_major, 0.025),
      p_ge1_major_p97_5 = qfun(p_ge1_major, 0.975),
      
      mean_q_import_mean  = mean(mean_q_import, na.rm = TRUE),
      mean_q_import_p2_5  = qfun(mean_q_import, 0.025),
      mean_q_import_p97_5 = qfun(mean_q_import, 0.975),
      
      mean_p_month_major_mean  = mean(mean_p_month_major, na.rm = TRUE),
      mean_p_month_major_p2_5  = qfun(mean_p_month_major, 0.025),
      mean_p_month_major_p97_5 = qfun(mean_p_month_major, 0.975),
      
      mean_p_est_mean  = mean(mean_p_est, na.rm = TRUE),
      mean_p_est_p2_5  = qfun(mean_p_est, 0.025),
      mean_p_est_p97_5 = qfun(mean_p_est, 0.975),
      .groups = "drop"
    )
  
  suffix <- sprintf("_%d_%d_rep%d", year_min, year_max, as.integer(n_rep))
  fp_monthly_csv <- file.path(out_dir, paste0("ctmc_spark_monthly",  suffix, ".csv"))
  fp_monthly_rds <- file.path(out_dir, paste0("ctmc_spark_monthly",  suffix, ".rds"))
  fp_yearly_csv  <- file.path(out_dir, paste0("ctmc_spark_yearly",   suffix, ".csv"))
  fp_yearly_rds  <- file.path(out_dir, paste0("ctmc_spark_yearly",   suffix, ".rds"))
  fp_horizon_csv <- file.path(out_dir, paste0("ctmc_spark_horizon",  suffix, ".csv"))
  fp_horizon_rds <- file.path(out_dir, paste0("ctmc_spark_horizon",  suffix, ".rds"))
  
  readr::write_csv(monthly_sum, fp_monthly_csv)
  saveRDS(monthly_sum, fp_monthly_rds)
  readr::write_csv(yearly_sum,  fp_yearly_csv)
  saveRDS(yearly_sum,  fp_yearly_rds)
  readr::write_csv(horizon_sum, fp_horizon_csv)
  saveRDS(horizon_sum, fp_horizon_rds)
  
  message("✅ Outputs written (summary over reps) to: ", out_dir)
  message("- monthly: ", fp_monthly_csv)
  message("- yearly:  ", fp_yearly_csv)
  message("- horizon: ", fp_horizon_csv)
  
  invisible(list(
    monthly = monthly_sum,
    yearly = yearly_sum,
    horizon = horizon_sum,
    monthly_reps = monthly_all,
    yearly_reps = yearly_all,
    horizon_reps = horizon_all
  ))
}

## ----------------------------------------------------------
## 8) Entrypoint
##    m = 1.0 (Ross-Macdonald standard)
## ----------------------------------------------------------
if (identical(environment(), globalenv())) {
  run_ctmc_spark(
    year_min = 2025,
    year_max = 2075,
    tau = 30L,
    infectious_period_days = 5,
    beta_vh = 0.3,
    beta_hv = 0.33,
    m = 1.0,         # Ross-Macdonald standard; sensitivity at {0.5, 0.8, 1.2, 2.0}
    include_import_in_q = FALSE,
    use_stochastic_EIP = TRUE,
    n_mc_eip = 2000L,
    n_rep = 1000L,
    seed = 123L,
    major_rule = list(type = "establishment"),
    out_dir = file.path(DIR_OUTPUT_SSP, "simulation")
  )
}

message("✓ ctmc_spark_monte_carlo.R loaded (SSP-aware + major_rule + p_local_major + numerical safeguards)")