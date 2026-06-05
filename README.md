# Dengue Establishment Risk in Turkey: CTMC Spark Model

## Project overview

This project models dengue establishment risk across five sentinel districts
in Turkey using a Continuous-Time Markov Chain (CTMC) spark model under
CMIP6 climate change scenarios (SSP1-2.6, SSP2-4.5, SSP5-8.5).

**Thesis**: "Aedes Kaynaklı Arbovirus Hassasiyeti"
**Author**: Okan Derin
**Institution**: Istanbul Medipol University

## Directory structure

```
r_project_tez/
├── R/
│   ├── 01_setup/      # Initialization, paths, packages, global options
│   ├── 02_data/       # Data import, cleaning, processing
│   ├── 03_models/     # CTMC spark, Monte Carlo, sensitivity analysis
│   ├── 04_results/    # Metrics derivation, Rmd reports
│   ├── 05_outputs/    # Figures, tables, maps
│   └── 06_simulation/ # (Future: cross-scenario comparisons)
├── data_raw/          # Original data (NOT tracked in git)
│   ├── climate/cmip6/{ssp126,ssp245,ssp585}/  # NetCDF files
│   ├── geography/     # Shapefiles
│   └── population/    # TurkStat, TUIK, tourism data
├── data_interim/      # Intermediate processed data
│   └── {ssp}/         # Scenario-specific
├── data_processed/    # Analysis-ready data
│   ├── {ssp}/         # Scenario-specific (climate, importation)
│   └── *.csv          # Shared (population, traits, species)
├── outputs/
│   ├── {ssp}/         # Per-scenario results
│   │   ├── model_results/
│   │   ├── sensitivity/
│   │   ├── figures/
│   │   └── tables/
│   └── cross_scenario/  # SSP comparison outputs
├── logs/              # Processing logs
├── renv.lock          # Package versions (reproducibility)
└── r_project_tez.Rproj
```

## Reproducibility

### Quick start

```r
# 1. Restore R package environment
renv::restore()

# 2. Run single scenario
Sys.setenv(SSP_SCENARIO = "ssp245")
source("R/run_all.R")

# 3. Run all scenarios
source("R/run_all.R")
```

### Requirements

- R >= 4.3.0
- RStudio (recommended)
- ~10 GB disk space for CMIP6 NetCDF data
- Required raw data files (see data_raw/README.md)

### Climate data

CMIP6 projections from EC-Earth3-CC model, downloaded via ESGF:
- SSP1-2.6, SSP2-4.5, SSP5-8.5
- Variables: tas, tasmax, tasmin, pr, hurs
- Resolution: Monthly, 2015–2100

## Key model parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| β_vh (vector→human) | 0.30 | Literature |
| β_hv (human→vector) | 0.33 | Literature |
| IP (infectious period) | 5 days | Literature |
| τ (epidemic threshold) | 30 | Model calibration |
| m (baseline mosquito:human) | 1.0 | Ross-Macdonald |

## License

This project is part of a master thesis. Please contact the author
before reusing data or code.

