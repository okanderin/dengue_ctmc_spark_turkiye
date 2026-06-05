## =========================================================
## Project : r_project_tez
## Script  : R/02_data/build_population_dataset.R
## Purpose : Sentinel district population (2015–2100)
##           Constant share method (baseline = 2024)
## =========================================================

## ---- Init (loads ROOT, packages, paths.R)
source("R/01_setup/init.R")

## ----------------------------
## User parameters
## ----------------------------
BASELINE_YEAR <- 2024
YEAR_MIN <- 2015
YEAR_MAX <- 2100
SCENARIO <- "main"   # main | low | high

## ----------------------------
## Logging
## ----------------------------
log_file <- file.path(DIR_LOGS, "build_population_dataset_log.txt")
sink(log_file, split = TRUE)
on.exit({ try(sink(), silent = TRUE) }, add = TRUE)

cat("=============================================\n")
cat("build_population_dataset.R started\n")
cat("Time          :", as.character(Sys.time()), "\n")
cat("Scenario      :", SCENARIO, "\n")
cat("Baseline year :", BASELINE_YEAR, "\n")
cat("Years         :", YEAR_MIN, "-", YEAR_MAX, "\n")
cat("=============================================\n\n")

## ----------------------------
## Helper functions
## ----------------------------
norm_tr <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- stringr::str_to_lower(x)
  stringi::stri_trans_general(x, "Latin-ASCII")
}

parse_num <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\s+", "")
  x <- stringr::str_replace_all(x, "\\.", "")
  x <- stringr::str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

## =========================================================
## 1) Sentinel baseline (2024)
## =========================================================
baseline_path <- file.path(DIR_POP_TURK, "sentinel_baseline_pop.xlsx")
stopifnot(file.exists(baseline_path))

baseline_raw <- readxl::read_excel(baseline_path)

baseline <- baseline_raw %>%
  dplyr::transmute(
    province_name = Province,
    district_name = District,
    year = BASELINE_YEAR,
    population = parse_num(pop_2024),
    prov_norm = norm_tr(Province),
    dist_norm = norm_tr(District)
  )

cat("--- Sentinel baseline ---\n")
print(baseline)

if (nrow(baseline) != 5) {
  stop("Sentinel baseline must contain exactly 5 rows.", call. = FALSE)
}

## =========================================================
## 2) GADM ADM2 lookup (district_id)
## =========================================================
shp_path <- file.path(DIR_SHP, "gadm41_TUR_2.shp")
stopifnot(file.exists(shp_path))

gadm2 <- sf::st_read(shp_path, quiet = TRUE)

gadm_lookup <- gadm2 %>%
  dplyr::mutate(
    NAME_1 = ifelse(NAME_1 == "Zinguldak", "Zonguldak", NAME_1),
    prov_norm = norm_tr(NAME_1),
    dist_norm = norm_tr(NAME_2),
    district_id = as.character(GID_2),
    province_name_gadm = NAME_1,
    district_name_gadm = NAME_2
  ) %>%
  dplyr::select(
    district_id,
    province_name_gadm,
    district_name_gadm,
    prov_norm,
    dist_norm
  ) %>%
  dplyr::distinct()

baseline_joined <- baseline %>%
  dplyr::left_join(gadm_lookup, by = c("prov_norm", "dist_norm"))

if (any(is.na(baseline_joined$district_id))) {
  print(baseline_joined)
  stop("Baseline → GADM join failed for some sentinel districts.", call. = FALSE)
}

baseline_final <- baseline_joined %>%
  dplyr::transmute(
    district_id,
    province_name = province_name_gadm,
    district_name = district_name_gadm,
    year,
    population
  )

## Save baseline
baseline_out <- file.path(DIR_PROCESSED, "population_sentinel_baseline_2024.csv")
write.csv(baseline_final, baseline_out, row.names = FALSE)
cat("Saved:", baseline_out, "\n\n")

## =========================================================
## 3) National population scenarios (TÜİK)
## =========================================================
scenario_path <- file.path(DIR_POP_PROJ, "pop_scenarios_tuik.xls")
stopifnot(file.exists(scenario_path))

scen_raw <- readxl::read_excel(scenario_path)

scenarios <- scen_raw %>%
  dplyr::transmute(
    year = as.integer(year),
    main = parse_num(main_scenario),
    low  = parse_num(low_scenario),
    high = parse_num(high_scenario)
  ) %>%
  dplyr::filter(year >= YEAR_MIN, year <= YEAR_MAX) %>%
  dplyr::arrange(year)

stopifnot(BASELINE_YEAR %in% scenarios$year)

turkey_total <- scenarios %>%
  dplyr::transmute(
    year,
    turkey_pop = dplyr::case_when(
      SCENARIO == "main" ~ main,
      SCENARIO == "low"  ~ low,
      SCENARIO == "high" ~ high
    )
  )

## =========================================================
## 4) Constant share projection (sentinel-only; growth-factor approach)
##   - Keep district shares constant relative to 2024 baseline
##   - Apply national growth factor over time
## =========================================================

# 1) National growth factor relative to baseline year
turkey_base <- turkey_total %>%
  dplyr::filter(year == BASELINE_YEAR) %>%
  dplyr::summarise(turkey_pop_base = dplyr::first(turkey_pop)) %>%
  dplyr::pull(turkey_pop_base)

if (length(turkey_base) == 0 || is.na(turkey_base)) {
  stop("Could not find BASELINE_YEAR in turkey_total for growth-factor.", call. = FALSE)
}

growth_tbl <- turkey_total %>%
  dplyr::mutate(growth_factor = turkey_pop / turkey_base)

# 2) Apply growth factor to each sentinel district's baseline population
population_yearly <- baseline_final %>%
  dplyr::transmute(
    district_id,
    province_name,
    district_name,
    pop_base = population
  ) %>%
  dplyr::cross_join(growth_tbl) %>%
  dplyr::mutate(
    population = pop_base * growth_factor,
    baseline_year = BASELINE_YEAR,
    scenario = SCENARIO
  ) %>%
  dplyr::select(
    district_id,
    province_name,
    district_name,
    year,
    population,
    scenario,
    baseline_year
  ) %>%
  dplyr::arrange(district_id, year)

# 3) Sanity check (baseline year) — should match exactly (within tiny tolerance)
check_2024 <- population_yearly %>%
  dplyr::filter(year == BASELINE_YEAR) %>%
  dplyr::left_join(
    baseline_final %>% dplyr::select(district_id, pop_base = population),
    by = "district_id"
  ) %>%
  dplyr::mutate(diff = population - pop_base)

if (max(abs(check_2024$diff), na.rm = TRUE) > 1e-6) {
  print(check_2024)
  stop("Baseline mismatch after growth-factor projection.", call. = FALSE)
}


## =========================================================
## 5) Save outputs
## =========================================================
out_csv <- file.path(DIR_PROCESSED, "population_sentinel_yearly_2015_2100.csv")
out_rds <- file.path(DIR_PROCESSED, "population_sentinel_yearly_2015_2100.rds")

write.csv(population_yearly, out_csv, row.names = FALSE)
saveRDS(population_yearly, out_rds)

cat("Saved:", out_csv, "\n")
cat("Saved:", out_rds, "\n")

cat("\n=============================================\n")
cat("build_population_dataset.R finished\n")
cat("Rows      :", nrow(population_yearly), "\n")
cat("Districts :", length(unique(population_yearly$district_id)), "\n")
cat("Years     :", min(population_yearly$year), "-", max(population_yearly$year), "\n")
cat("=============================================\n")
