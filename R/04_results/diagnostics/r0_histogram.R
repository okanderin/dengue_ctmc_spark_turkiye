library(readr); library(dplyr); library(ggplot2)

# --- Veriyi oku ---
ssp_map <- c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5", ssp585 = "SSP5-8.5")

r0 <- bind_rows(lapply(names(ssp_map), function(s) {
  read_csv(file.path("outputs", s, "simulation",
                     "ctmc_spark_monthly_2025_2075_rep1000.csv"), show_col_types = FALSE) |>
    transmute(R0 = lambda_local_i1_mean / 0.2,
              Senaryo = ssp_map[[s]])
})) |> filter(R0 > 0)

# --- Panel etiketlerine R0>1 sayısını göm (dinamik) ---
etiket <- r0 |>
  group_by(Senaryo) |>
  summarise(n_sup = sum(R0 > 1), .groups = "drop") |>
  mutate(lab = sprintf("%s (R\u2080>1: %d ay, %%%.1f)", Senaryo, n_sup, 100 * n_sup / 3060))

r0 <- r0 |>
  left_join(etiket, by = "Senaryo") |>
  mutate(panel = factor(lab, levels = etiket$lab))

# --- Grafik ---
ggplot(r0, aes(R0)) +
  geom_histogram(bins = 60, fill = "grey35") +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .7, color = "red") +
  annotate("text", x = 1.05, y = Inf, label = "R\u2080 = 1", hjust = 0, vjust = 2,
           color = "red", size = 3) +
  facet_wrap(~ panel) +
  labs(title = "Yerel üreme sayısının (R\u2080) dağılımı — aktif aylar",
       subtitle = "Kesikli çizgi: süperkritik eşik (R\u2080 = 1). Toplam 3.060 ilçe-ay (5 ilçe × 51 yıl × 12 ay).",
       x = expression(R[0] ~ "(yerel)"), y = "Sıklık") +
  theme_minimal(base_size = 11)

ggsave("outputs/cross_scenario/fig_R0_dagilim.png", width = 11, height = 4, dpi = 300)
cat("✓ Kaydedildi: outputs/cross_scenario/fig_R0_dagilim.png\n")
print(etiket)