# ==========================================================
# packages.R
# Project : r_project_tez
# Purpose : Centralized package loader (NO auto-install)
# Author  : Okan Derin
# ==========================================================

## ----------------------------------------------------------
## 1) Package groups (logical blocks)
## ----------------------------------------------------------

pkgs <- c(
  
  # ---- Core data wrangling ----
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "lubridate",
  "tibble",
  "purrr",
  "tidyverse",
  "gridExtra",
  
  # ---- Text / encoding robustness (Turkish chars etc.) ----
  "stringi",
  
  # ---- Spatial & raster (climate, districts) ----
  "sf",
  "terra",
  "exactextractr",
  
  # ---- Project structure & reproducibility ----
  "here",
  
  # ---- Diagnostics / utilities ----
  "rlang",
  
  # CTMC / stochastic phase
  "expm", "Matrix", 
  
  # Visualization
  "ggplot2", "patchwork", "igraph", "tmap", "grid",
  
  # Tables / reporting
  "gtsummary", "flextable"
  
)

pkgs <- unique(pkgs)

## ----------------------------------------------------------
## 2) Check missing packages
## ----------------------------------------------------------

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  stop(
    "\n=================================================\n",
    "Missing required R packages detected:\n\n  - ",
    paste(missing, collapse = "\n  - "),
    "\n\nPlease install them ONCE, then re-run the project.\n\n",
    "Base R:\n",
    "  install.packages(c(",
    paste(sprintf('\"%s\"', missing), collapse = ", "),
    "))\n\n",
    "renv (recommended for thesis reproducibility):\n",
    "  renv::install(c(",
    paste(sprintf('\"%s\"', missing), collapse = ", "),
    "))\n",
    "=================================================\n",
    call. = FALSE
  )
}

## ----------------------------------------------------------
## 3) Load packages quietly
## ----------------------------------------------------------

invisible(
  lapply(pkgs, function(p) {
    suppressPackageStartupMessages(
      library(p, character.only = TRUE)
    )
  })
)

## ----------------------------------------------------------
## 4) Minimal startup confirmation (optional)
## ----------------------------------------------------------

message("✔ All required packages loaded successfully.")
