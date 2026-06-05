# ==========================================================
# DIAGNOSTIC ANALYSIS
# Transmission mechanism validation for CTMC spark model
# ==========================================================

source("R/01_setup/init.R")

library(dplyr)
library(ggplot2)
library(readr)

# ----------------------------------------------------------
# 1 Load monthly model results
# ----------------------------------------------------------

fp_monthly <- file.path("outputs/model_results/ctmc_spark_monthly_2025_2075.csv")

monthly <- read_csv(fp_monthly, show_col_types = FALSE)

# ----------------------------------------------------------
# 2 Validate required columns
# ----------------------------------------------------------

required_monthly <- c(
  "district_id",
  "year",
  "month",
  "temp_c",
  "lambda_import",
  "lambda_local_i1"
)

missing_monthly <- setdiff(required_monthly, names(monthly))

if (length(missing_monthly) > 0) {
  stop("Missing required columns in monthly: ", paste(missing_monthly, collapse = ", "))
}

# ----------------------------------------------------------
# 3 Create R0_local if missing
# ----------------------------------------------------------

if (!"R0_local" %in% names(monthly)) {
  
  message("R0_local column missing — computing from lambda_local_i1")
  
  gamma_h <- 1/5   # human recovery rate
  
  monthly <- monthly %>%
    mutate(
      R0_local = lambda_local_i1 / gamma_h
    )
  
}

# ----------------------------------------------------------
# 4 Distribution diagnostics
# ----------------------------------------------------------

dist_summary <- monthly %>%
  summarise(
    min_R0 = min(R0_local, na.rm = TRUE),
    median_R0 = median(R0_local, na.rm = TRUE),
    mean_R0 = mean(R0_local, na.rm = TRUE),
    max_R0 = max(R0_local, na.rm = TRUE)
  )

print(dist_summary)

# ----------------------------------------------------------
# 5 Plot 1 — R0 distribution
# ----------------------------------------------------------

p1 <- ggplot(monthly, aes(x = R0_local)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "Distribution of Local Reproduction Number",
    x = "R0_local",
    y = "Frequency"
  ) +
  theme_minimal()

# ----------------------------------------------------------
# 6 Plot 2 — Temperature vs R0
# ----------------------------------------------------------

p2 <- ggplot(monthly, aes(temp_c, R0_local)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Temperature vs Local R0",
    x = "Temperature (°C)",
    y = "R0_local"
  ) +
  theme_minimal()

# ----------------------------------------------------------
# 7 Plot 3 — Temperature vs lambda_local
# ----------------------------------------------------------

p3 <- ggplot(monthly, aes(temp_c, lambda_local_i1)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Temperature vs Local Transmission Rate",
    x = "Temperature (°C)",
    y = "lambda_local"
  ) +
  theme_minimal()

# ----------------------------------------------------------
# 8 Plot 4 — Importation vs establishment
# ----------------------------------------------------------

if ("p_establishment" %in% names(monthly)) {
  
  p4 <- ggplot(monthly, aes(lambda_import, p_establishment)) +
    geom_point(alpha = 0.25) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      title = "Importation Pressure vs Establishment Probability",
      x = "lambda_import",
      y = "p_establishment"
    ) +
    theme_minimal()
  
}

# ----------------------------------------------------------
# 9 Save results
# ----------------------------------------------------------

dir.create("outputs/diagnostics", showWarnings = FALSE)

ggsave("outputs/diagnostics/R0_distribution.png", p1, width = 6, height = 4)
ggsave("outputs/diagnostics/temperature_vs_R0.png", p2, width = 6, height = 4)
ggsave("outputs/diagnostics/temperature_vs_lambda_local.png", p3, width = 6, height = 4)

if (exists("p4")) {
  ggsave("outputs/diagnostics/import_vs_establishment.png", p4, width = 6, height = 4)
}

write_csv(dist_summary, "outputs/diagnostics/R0_summary.csv")

message("Diagnostics completed successfully")