# ==========================================================
# Seasonality + warming diagnostics for CTMC spark outputs
# ==========================================================

source("R/01_setup/init.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(forcats)
  library(scales)
})

# ----------------------------------------------------------
# 0) Read monthly results
# ----------------------------------------------------------
fp_monthly <- file.path("outputs", "model_results", "ctmc_spark_monthly_2025_2075.rds")
monthly <- readRDS(fp_monthly)

dir_out <- file.path("outputs", "diagnostics", "ctmc_spark")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Eğer R0_local yoksa üret
if (!"R0_local" %in% names(monthly)) {
  gamma_h <- 1 / 5
  monthly <- monthly %>%
    mutate(R0_local = lambda_local_i1 / gamma_h)
}

# ----------------------------------------------------------
# 1) Seasonality table
# ----------------------------------------------------------
month_labs_tr <- c("Oca", "Şub", "Mar", "Nis", "May", "Haz",
                   "Tem", "Ağu", "Eyl", "Eki", "Kas", "Ara")

tab_seasonality <- monthly %>%
  group_by(month) %>%
  summarise(
    n = n(),
    mean_temp = mean(temp_c, na.rm = TRUE),
    mean_rh = mean(rh, na.rm = TRUE),
    mean_lambda_local = mean(lambda_local_i1, na.rm = TRUE),
    median_lambda_local = median(lambda_local_i1, na.rm = TRUE),
    max_lambda_local = max(lambda_local_i1, na.rm = TRUE),
    mean_R0 = mean(R0_local, na.rm = TRUE),
    median_R0 = median(R0_local, na.rm = TRUE),
    max_R0 = max(R0_local, na.rm = TRUE),
    mean_p_est = if ("p_establishment" %in% names(.)) mean(p_establishment, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    month_f = factor(month, levels = 1:12, labels = month_labs_tr)
  )

write_csv(tab_seasonality, file.path(dir_out, "seasonality_summary.csv"))

# ----------------------------------------------------------
# 2) Plot: seasonality of R0
# ----------------------------------------------------------
p_season_r0 <- tab_seasonality %>%
  ggplot(aes(x = month_f, y = mean_R0, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.8) +
  labs(
    title = "Türkiye'de dang uygunluğunun mevsimselliği",
    subtitle = "Aylık ortalama yerel çoğalma sayısı (R0_local)",
    x = "Ay",
    y = "Ortalama R0_local"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(dir_out, "fig_seasonality_R0_local.png"),
  p_season_r0, width = 9, height = 5, dpi = 320
)

# ----------------------------------------------------------
# 3) Plot: seasonality of lambda_local
# ----------------------------------------------------------
p_season_lambda <- tab_seasonality %>%
  ggplot(aes(x = month_f, y = mean_lambda_local, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Türkiye'de dang uygunluğunun mevsimselliği",
    subtitle = "Aylık ortalama yerel transmisyon hızı (lambda_local_i1)",
    x = "Ay",
    y = "Ortalama lambda_local_i1"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(dir_out, "fig_seasonality_lambda_local.png"),
  p_season_lambda, width = 9, height = 5, dpi = 320
)

# ----------------------------------------------------------
# 4) Plot: seasonality of establishment probability
# ----------------------------------------------------------
if ("p_establishment" %in% names(monthly)) {
  p_season_pest <- tab_seasonality %>%
    ggplot(aes(x = month_f, y = mean_p_est, group = 1)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = "Türkiye'de dang uygunluğunun mevsimselliği",
      subtitle = "Aylık ortalama establishment olasılığı",
      x = "Ay",
      y = "Ortalama p_establishment"
    ) +
    theme_minimal(base_size = 13)
  
  ggsave(
    file.path(dir_out, "fig_seasonality_p_establishment.png"),
    p_season_pest, width = 9, height = 5, dpi = 320
  )
}

# ----------------------------------------------------------
# 5) District-specific seasonality
# ----------------------------------------------------------
tab_seasonality_district <- monthly %>%
  group_by(district_id, month) %>%
  summarise(
    mean_temp = mean(temp_c, na.rm = TRUE),
    mean_lambda_local = mean(lambda_local_i1, na.rm = TRUE),
    mean_R0 = mean(R0_local, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(month_f = factor(month, levels = 1:12, labels = month_labs_tr))

write_csv(tab_seasonality_district, file.path(dir_out, "seasonality_by_district.csv"))

p_season_r0_district <- tab_seasonality_district %>%
  ggplot(aes(x = month_f, y = mean_R0, group = district_id, color = district_id)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.8) +
  labs(
    title = "İlçe bazında dang uygunluğunun mevsimselliği",
    subtitle = "Aylık ortalama R0_local",
    x = "Ay",
    y = "Ortalama R0_local",
    color = "İlçe"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(dir_out, "fig_seasonality_R0_local_by_district.png"),
  p_season_r0_district, width = 10, height = 6, dpi = 320
)

# ----------------------------------------------------------
# 6) Climate warming scenario: direct temperature-shift probe
#    Not a full rerun. Exploratory perturbation analysis.
# ----------------------------------------------------------

warming_grid <- c(0, 1, 2, 3, 4)

# Bu yaklaşım yalnızca mevcut temp-R0 ilişkisini kaydırarak yaklaşık bir gösterim üretir.
# Tam mekanistik senaryo için modelin tekrar koşulması gerekir.
#
# Burada observed temp -> observed R0 ilişkisini loess ile fit ediyoruz,
# sonra temp + delta için beklenen R0'ı tahmin ediyoruz.

loess_r0 <- loess(R0_local ~ temp_c, data = monthly, span = 0.5)
loess_lambda <- loess(lambda_local_i1 ~ temp_c, data = monthly, span = 0.5)

warming_df <- expand.grid(
  delta_temp = warming_grid,
  row_id = seq_len(nrow(monthly))
) %>%
  as_tibble() %>%
  left_join(
    monthly %>%
      mutate(row_id = row_number()) %>%
      select(row_id, district_id, year, month, temp_c, rh, lambda_local_i1, R0_local),
    by = "row_id"
  ) %>%
  mutate(
    temp_shifted = temp_c + delta_temp,
    R0_shifted = pmax(0, predict(loess_r0, newdata = data.frame(temp_c = temp_shifted))),
    lambda_shifted = pmax(0, predict(loess_lambda, newdata = data.frame(temp_c = temp_shifted)))
  )

tab_warming_summary <- warming_df %>%
  group_by(delta_temp) %>%
  summarise(
    mean_R0_shifted = mean(R0_shifted, na.rm = TRUE),
    median_R0_shifted = median(R0_shifted, na.rm = TRUE),
    max_R0_shifted = max(R0_shifted, na.rm = TRUE),
    prop_R0_gt1_shifted = mean(R0_shifted > 1, na.rm = TRUE),
    mean_lambda_shifted = mean(lambda_shifted, na.rm = TRUE),
    max_lambda_shifted = max(lambda_shifted, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(tab_warming_summary, file.path(dir_out, "warming_shift_summary.csv"))

# ----------------------------------------------------------
# 7) Plot: warming vs mean R0
# ----------------------------------------------------------
p_warming_r0 <- tab_warming_summary %>%
  ggplot(aes(x = delta_temp, y = mean_R0_shifted)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.8) +
  scale_x_continuous(breaks = warming_grid) +
  labs(
    title = "İklim ısınması ile yerel dang uygunluğundaki kayma",
    subtitle = "Sıcaklık +Δ°C kaydırma altında ortalama R0_local",
    x = expression(Delta*" Sıcaklık (°C)"),
    y = "Ortalama R0_local"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(dir_out, "fig_warming_shift_mean_R0.png"),
  p_warming_r0, width = 8, height = 5, dpi = 320
)

# ----------------------------------------------------------
# 8) Plot: warming vs max R0
# ----------------------------------------------------------
p_warming_max_r0 <- tab_warming_summary %>%
  ggplot(aes(x = delta_temp, y = max_R0_shifted)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.8) +
  scale_x_continuous(breaks = warming_grid) +
  labs(
    title = "İklim ısınması ile pik yerel dang uygunluğundaki kayma",
    subtitle = "Sıcaklık +Δ°C kaydırma altında maksimum R0_local",
    x = expression(Delta*" Sıcaklık (°C)"),
    y = "Maksimum R0_local"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(dir_out, "fig_warming_shift_max_R0.png"),
  p_warming_max_r0, width = 8, height = 5, dpi = 320
)

# ----------------------------------------------------------
# 9) Plot: warming vs proportion above threshold
# ----------------------------------------------------------
p_warming_prop <- tab_warming_summary %>%
  ggplot(aes(x = delta_temp, y = prop_R0_gt1_shifted)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = warming_grid) +
  labs(
    title = "İklim ısınması altında eşik üstü yerel uygunluk",
    subtitle = "Sıcaklık +Δ°C kaydırma altında R0_local > 1 oranı",
    x = expression(Delta*" Sıcaklık (°C)"),
    y = "R0_local > 1 oranı"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(dir_out, "fig_warming_shift_prop_R0_gt1.png"),
  p_warming_prop, width = 8, height = 5, dpi = 320
)

# ----------------------------------------------------------
# 10) Thesis note
# ----------------------------------------------------------
note_lines <- c(
  "Seasonality and warming-shift note",
  "----------------------------------",
  "Seasonality figures are descriptive summaries of observed model outputs.",
  "Warming-shift figures are exploratory perturbation analyses based on the empirical temperature-R0 and temperature-lambda relationships.",
  "They do NOT replace a full mechanistic rerun under IPCC scenario-specific climate inputs."
)

writeLines(note_lines, con = file.path(dir_out, "seasonality_warming_note.txt"))

message("Seasonality and warming-shift diagnostics completed.")
message("Outputs written to: ", normalizePath(dir_out))