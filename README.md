# Dengue establishment risk in Türkiye under CMIP6 scenarios — reproducibility package

R pipeline for a **continuous-time Markov chain (CTMC) spark-phase birth–death model** that estimates the probability an imported dengue infection establishes a local transmission chain (reaching an operational threshold of τ = 30 concurrent infectious humans) across five climatically distinct sentinel districts in Türkiye, under CMIP6 **SSP1-2.6, SSP2-4.5 and SSP5-8.5** scenarios over **2025–2075**.

This repository accompanies the manuscript *"A stochastic model of dengue establishment following importation into Türkiye under CMIP6 climate scenarios"* and its Supporting Information (S1–S5 Files, S1 Data).

Sentinel districts: **Kartal/İstanbul** (*Ae. albopictus*), **Fethiye/Muğla** (*Ae. albopictus*), **Hopa/Artvin** (*Ae. aegypti*), **Zonguldak Merkez** (*Ae. aegypti* frontier scenario), **Eğirdir/Isparta** (*Ae. albopictus*).

---

## Repository layout

```
.
├── R/                     # analysis pipeline (numbered by stage)
│   ├── 01_setup/          # init.R, packages.R, paths.R, global_options.R
│   ├── 02_data/           # climate bias-correction, importation pressure, trait params, travel weights
│   ├── 03_models/         # CTMC spark model, nested Monte Carlo, parameter functions, sensitivity
│   ├── 04_results/        # per-SSP generation, figures, regression, verification & check scripts
│   ├── run_all.R          # end-to-end driver
│   └── data_build_once_run.R
├── data_processed/        # non-proprietary processed inputs (bias-corrected climate summaries,
│                          #   trait_params_*.csv, travel_weights_*, seasonality, sentinel_species, GBD inputs)
├── data_raw/              # small public inputs + README documenting CDS queries for raw NetCDF (NOT redistributed)
├── outputs/               # district–month–year model outputs, tables, validation, figures
│   ├── sspXXX/simulation/ # ctmc_spark_{monthly,yearly,horizon}_2025_2075_rep1000.rds  (per scenario)
│   ├── sspXXX/validation/ # Stage 1 / Stage 2 (Jensen) verification
│   └── tables/            # cross-scenario tables (regression, PRCC, T_opt, Poisson–Bernoulli, k_vpd)
├── shiny_dengue_app/      # companion Shiny app (projection browser, What-If, arbitrary-location risk)
├── renv.lock              # locked package environment (see "Reproducing")
├── CITATION.cff
├── LICENSE
└── README.md
```

## Model, in brief

- **Spark-phase CTMC** on the number of infectious humans I ∈ {0, …, τ}, τ = 30; establishment probability from the finite-threshold gambler's-ruin formula (S1 File, Eq. 1).
- **Local transmission rate** λ(i=1)(T,RH) via a reduced Ross–Macdonald vectorial-capacity construction; R₀ = λ(i=1)/γ, γ = 0.20 day⁻¹ (S1 File, Eq. 2).
- **Temperature-/humidity-dependent vector biology**: Brière biting rate and EIP development, quadratic lifespan, and a VPD-based mortality multiplier μ_v = (1/lf(T))·exp(k_vpd·(VPD − VPD_ref)), VPD_ref = 1.0 kPa, k_vpd = 0.5 kPa⁻¹ (S1/S2 Files, Eq. 3).
- **EIP heterogeneity + Jensen correction**: individual EIP ~ LogNormal at fixed T; nested Monte Carlo (inner n_mc = 2,000) estimates E[P_est] rather than P_est(E[EIP]) (S1 File, Eq. 4).
- **Importation pressure**: province-level foreign-arrival volume × GBD 2023 travel-weighted source-country incidence × seasonality × viremia-at-arrival, as a non-homogeneous Poisson process (S2 File, Eq. 5a–5g).
- **Outer replication** n_rep = 1,000 is seed-controlled reproducibility, *not* parametric uncertainty; parametric uncertainty is assessed separately by LHS–PRCC and OAT (S4 File).

## Pipeline order

```r
# From the repository root, with the R project open:
source("R/01_setup/init.R")          # loads packages.R + paths.R + options
source("R/data_build_once_run.R")    # 02_data: bias correction, importation pressure, trait params, ...
source("R/run_all.R")                # 03_models + 04_results: simulate all SSPs, tables, figures
```

Key stage scripts: `R/03_models/ctmc_spark.R`, `ctmc_spark_monte_carlo.R`, `parameter_functions.R`; results in `R/04_results/01_generate_ssp_outputs.R`, `02_all_ssp_generate.R`, `bulgular.Rmd`.

## Verification and sensitivity

- **Two-stage verification** (`R/run_all.R` → `validation_two_stage*.R`, `R/04_results/mc_validation_vs_analytic.R`): Stage 1 reproduces the analytic finite-threshold solution to machine precision (σ_EIP = 0) and against an independent brute-force CTMC; Stage 2 quantifies the Jensen bias. Outputs in `outputs/sspXXX/validation/`.
- **Poisson–binomial check** (`R/04_results/ek_a4_poisson_bernoulli_check.R`): regenerates Appendix Table A.4 — the exact heterogeneous-Bernoulli complement product vs the Poisson cumulative-hazard approximation — from the monthly outputs. Result: `outputs/tables/tbl_ek_a4_poisson_bernoulli.csv`.
- **k_vpd sensitivity** (`R/04_results/kvpd_sensitivity_check.R`): recomputes P_est from each district's actual T/RH trajectory at k_vpd = 0.3, 0.5, 0.8 kPa⁻¹ under all three scenarios, and reports the effect on P_horizon and the district ranking. Results: `outputs/tables/tbl_kvpd_sensitivity.csv`, `tbl_kvpd_reference_check.csv`.
- **LHS–PRCC / OAT** (`R/03_models/sensitivity_ctmc_mc.R`, `LHS_PRCC_EN_ref_1-1.R`, `OAT_*`): global and one-at-a-time sensitivity.

## Reproducing the environment

Package versions are pinned with [`renv`](https://rstudio.github.io/renv/):

```r
install.packages("renv")
renv::restore()      # installs the exact package versions from renv.lock
```

R ≥ 4.3 is assumed. A full run takes roughly 45–60 hours on an Intel i7-6700HQ / 16 GB RAM.

## Data availability

- **Processed inputs** used by the pipeline are in `data_processed/` (CC BY 4.0).
- **Raw CMIP6 (CNRM-CM6-1-HR) and ERA5-Land NetCDF** are **not redistributed** because of size and licensing; `data_raw/README_data.md` documents the exact Copernicus Climate Data Store queries needed to retrieve them.
- **GBD 2023** dengue incidence inputs are included under `data_raw/population/GBD_2023_DATA/` with the IHME citation.

## Citation

See `CITATION.cff`. Please cite both this software and the associated PLOS NTD article.

## Licence

Code: MIT. Processed data and outputs: CC BY 4.0. Raw climate data: terms of the original providers. See `LICENSE`.
