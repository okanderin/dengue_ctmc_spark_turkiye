#source("R/01_setup/init.R")
# # -------------------------------------------------------------------
# init.R
# Project initializer for r_project_tez
# Guarantees that DIR_* objects are visible in analysis scripts
# -------------------------------------------------------------------


.source_safe <- function(path, envir) {
  if (!file.exists(path)) {
    stop("Missing init dependency: ", path, call. = FALSE)
  }
  source(path, local = envir)
}

# --- Require here
if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. Install once: install.packages('here')",
    call. = FALSE
  )
}

# --- Determine ROOT
ROOT <- here::here()

if (basename(ROOT) != "r_project_tez") {
  candidate <- file.path(ROOT, "r_project_tez")
  if (dir.exists(candidate)) {
    ROOT <- candidate
  } else {
    stop(
      "Project root could not be validated as 'r_project_tez'.\n",
      "Detected ROOT: ", ROOT,
      call. = FALSE
    )
  }
}

SETUP_DIR <- file.path(ROOT, "R", "01_setup")
if (!dir.exists(SETUP_DIR)) {
  stop("Setup directory not found: ", SETUP_DIR, call. = FALSE)
}

# 🔑 CRITICAL FIX:
# Always source setup scripts into the GLOBAL environment
TARGET_ENV <- globalenv()

.source_safe(file.path(SETUP_DIR, "packages.R"),       envir = TARGET_ENV)
.source_safe(file.path(SETUP_DIR, "paths.R"),          envir = TARGET_ENV)
.source_safe(file.path(SETUP_DIR, "global_options.R"), envir = TARGET_ENV)
