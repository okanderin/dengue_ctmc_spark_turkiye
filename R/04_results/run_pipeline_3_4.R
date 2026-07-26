## =========================================================
## run_pipeline_3_4.R
## Parallel wrapper for pipeline steps 3-4:
##   STEP 3  CTMC Monte Carlo per SSP        (parallel, isolated workers)
##   STEP 4a per-SSP generation (figures/tables from CTMC outputs)
##   STEP 4b cross-scenario figures          (needs all 3 SSPs)
##
## Windows-safe: uses future multisession (NOT mclapply; no fork on Windows).
## Each SSP runs in its own R session, so SSP_SCENARIO env var and the
## per-session EIP cache are fully isolated.
##
## Determinism: run_ctmc_spark seeds EVERY cell via
##   set.seed(base_seed + hash(district_id, year, month, rep_id))
## so results are independent of execution order. Parallel == sequential,
## bit for bit. (future.seed = TRUE only quiets the RNG warning; it does
## not change the per-cell seeded outputs.)
##
## Prereq: STEP 1-2 already produced travel_weights_zone.rds and the
##   zone-based importation_pressure_*.rds for each SSP.
##
## Usage (from project root):
##   source("R/run_pipeline_3_4.R")
## Optional cores:
##   Sys.setenv(N_WORKERS = 3); source("R/run_pipeline_3_4.R")
## =========================================================

suppressPackageStartupMessages({
  if (!requireNamespace("rprojroot",   quietly = TRUE)) install.packages("rprojroot")
  if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
  library(future.apply)
})

PROJ <- tryCatch(
  rprojroot::find_root(rprojroot::has_file("r_project_tez.Rproj")),
  error = function(e) getwd()
)
SSPS      <- c("ssp126", "ssp245", "ssp585")
CTMC_FILE <- file.path("R", "03_models", "ctmc_spark_monte_carlo.R")
GEN_FILE  <- file.path("R", "04_results", "01_generate_ssp_outputs.R")
XSCN_FILE <- file.path("R", "04_results", "02_figures_cross_scenario.R")
N_WORKERS <- as.integer(Sys.getenv("N_WORKERS",
                                   unset = as.character(min(3L, length(SSPS)))))

setwd(PROJ)
cat("=============================================\n")
cat("Pipeline 3-4 | PROJ =", PROJ, "\n")
cat("Workers =", N_WORKERS, " | SSPs =", paste(SSPS, collapse = ", "), "\n")
cat("=============================================\n")

## ---------------------------------------------------------
## STEP 3 — CTMC per SSP (parallel)
##   Sourcing ctmc_spark_monte_carlo.R at top level fires its
##   entrypoint (run_ctmc_spark) for the active SSP_SCENARIO.
## ---------------------------------------------------------
plan(multisession, workers = N_WORKERS)
on.exit(plan(sequential), add = TRUE)

t0 <- Sys.time()
step3 <- future_lapply(SSPS, function(s) {
  setwd(PROJ)                          # workers start fresh; pin the wd
  Sys.setenv(SSP_SCENARIO = s)
  ok <- tryCatch({
    source(CTMC_FILE, local = FALSE)   # local=FALSE -> runs in globalenv -> entrypoint fires
    "OK"
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  paste0(s, ": ", ok)
}, future.seed = TRUE)

plan(sequential)                       # release workers before STEP 4
cat("\n--- STEP 3 sonuclari ---\n"); print(unlist(step3))
cat("STEP 3 (CTMC x", length(SSPS), ") sure: ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " dk\n", sep = "")

if (any(grepl("ERROR", unlist(step3)))) {
  stop("STEP 3'te en az bir SSP basarisiz — STEP 4 atlandi.", call. = FALSE)
}

## ---------------------------------------------------------
## STEP 4a — per-SSP generation (reads CTMC outputs; fast, sequential)
## ---------------------------------------------------------
cat("\n--- STEP 4a: per-SSP generation ---\n")
for (s in SSPS) {
  setwd(PROJ)
  Sys.setenv(SSP_SCENARIO = s)
  cat(">> 01_generate_ssp_outputs.R  [", s, "]\n", sep = "")
  source(GEN_FILE, local = FALSE)
}

## ---------------------------------------------------------
## STEP 4b — cross-scenario figures (needs all 3 SSPs)
## ---------------------------------------------------------
if (file.exists(XSCN_FILE)) {
  cat("\n--- STEP 4b: cross-scenario figures ---\n")
  setwd(PROJ)
  source(XSCN_FILE, local = FALSE)
}

cat("\n=============================================\n")
cat("Pipeline 3-4 TAMAMLANDI.\n")
cat("Sonraki: 6.10_pest_regression.R (T_opt) + 01_generate P_ufuk tablosu.\n")
cat("=============================================\n")

## ---------------------------------------------------------
## Daha fazla cekirdek gerekiyorsa (opsiyonel):
##   SSP-seviyesi paralellik en fazla 3x hizlandirir. Cok cekirdekli
##   makinede rep-seviyesi paralellik daha fazla kazandirir; bunun icin
##   ctmc_spark_monte_carlo.R icindeki
##     reps <- lapply(rep_ids, compute_one_rep)
##   satirini sunla degistirin (Windows-uyumlu):
##     library(future.apply); plan(multisession)
##     reps <- future_lapply(rep_ids, compute_one_rep, future.seed = TRUE)
##   Per-hucre seed'ler yuzunden sonuc yine deterministiktir.
## ---------------------------------------------------------
