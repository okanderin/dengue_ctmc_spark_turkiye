# ==========================================================
# R/04_results/fig_decade_progression_3ssp.R
# Onyıllık otokton salgın riski progresyonu — ÜÇ SSP senaryosu (facet)
# - Log10 y-ekseni (değerler 10^-2 ... 10^-13, ~11 mertebe)
# - 3 panel yan yana, ortak y-ekseni -> senaryolar arası fark görünür
# - Çizgi + nokta; monotonik olmayan onyıllık deseni gösterir
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# ---- Veri (Tablo 6.6.1) -----------------------------------
df <- tibble::tribble(
  ~ssp,        ~ilce,       ~`2020`,   ~`2030`,   ~`2040`,   ~`2050`,   ~`2060`,   ~`2070`,
  "SSP1-2.6",  "Kartal",    3.10e-3,   5.20e-3,   4.68e-3,   8.05e-3,   9.73e-3,   7.21e-3,
  "SSP1-2.6",  "Fethiye",   1.22e-7,   1.47e-7,   1.04e-7,   4.85e-7,   2.44e-7,   9.15e-8,
  "SSP1-2.6",  "Zonguldak", 3.90e-9,   5.94e-8,   4.57e-7,   1.50e-7,   5.60e-7,   8.44e-7,
  "SSP1-2.6",  "Hopa",      3.40e-8,   1.03e-7,   2.59e-6,   7.17e-7,   3.69e-6,   1.16e-5,
  "SSP1-2.6",  "Eğirdir",   6.38e-13,  3.02e-12,  4.79e-13,  5.99e-12,  2.98e-11,  2.42e-12,
  "SSP2-4.5",  "Kartal",    3.98e-3,   3.50e-3,   7.49e-3,   9.99e-3,   1.21e-2,   1.56e-2,
  "SSP2-4.5",  "Fethiye",   2.83e-6,   3.76e-7,   1.97e-6,   4.45e-8,   1.54e-6,   3.81e-5,
  "SSP2-4.5",  "Zonguldak", 3.13e-7,   1.46e-8,   1.01e-6,   8.69e-7,   3.88e-6,   4.35e-6,
  "SSP2-4.5",  "Hopa",      9.87e-6,   2.61e-7,   9.71e-6,   5.42e-6,   4.29e-5,   6.65e-5,
  "SSP2-4.5",  "Eğirdir",   6.49e-12,  1.24e-12,  3.88e-11,  2.05e-12,  2.82e-11,  9.07e-10,
  "SSP5-8.5",  "Kartal",    6.44e-3,   6.54e-3,   9.67e-3,   1.32e-2,   2.18e-2,   1.85e-2,
  "SSP5-8.5",  "Fethiye",   1.08e-7,   7.34e-7,   2.15e-7,   3.86e-7,   1.05e-6,   1.17e-6,
  "SSP5-8.5",  "Zonguldak", 5.63e-7,   1.03e-7,   8.21e-7,   3.79e-6,   8.56e-6,   1.04e-5,
  "SSP5-8.5",  "Hopa",      9.13e-6,   6.68e-6,   1.05e-5,   3.83e-5,   8.58e-5,   1.37e-4,
  "SSP5-8.5",  "Eğirdir",   1.16e-12,  1.79e-11,  8.05e-12,  3.27e-11,  7.33e-11,  2.02e-11
)

df_long <- df %>%
  pivot_longer(c(`2020`,`2030`,`2040`,`2050`,`2060`,`2070`),
               names_to = "onyil", values_to = "risk") %>%
  mutate(
    onyil = paste0(onyil, "'ler"),
    onyil = factor(onyil, levels = c("2020'ler","2030'lar","2040'lar",
                                     "2050'ler","2060'lar","2070'ler")),
    ssp   = factor(ssp, levels = c("SSP1-2.6","SSP2-4.5","SSP5-8.5")),
    ilce  = factor(ilce, levels = c("Kartal","Hopa","Fethiye","Zonguldak","Eğirdir"))
  )

# ---- Renk / şekil (Cove kategorik) ------------------------
pal <- c("Kartal"="#2a78d6","Hopa"="#eda100","Fethiye"="#eb6834",
         "Zonguldak"="#1baf7a","Eğirdir"="#4a3aa7")
shp <- c("Kartal"=16,"Hopa"=18,"Fethiye"=15,"Zonguldak"=17,"Eğirdir"=25)

# ---- Grafik ----------------------------------------------
p <- ggplot(df_long, aes(onyil, risk, color = ilce, group = ilce, shape = ilce)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2, fill = "white", stroke = 0.7) +
  facet_wrap(~ ssp, nrow = 1) +
  scale_y_log10(
    breaks = 10^seq(-13, -2, 1),
    labels = trans_format("log10", math_format(10^.x)),
    minor_breaks = as.vector(outer(1:9, 10^seq(-14, -2, 1)))
  ) +
  scale_color_manual(values = pal, name = "İlçe") +
  scale_shape_manual(values = shp, name = "İlçe") +
  labs(
    x = "Onyıl",
    y = expression("Onyıllık otokton salgın riski"~(log[10]~"ölçek)")),
    title = "Onyıllık Risk Progresyonu — Üç SSP Senaryosu"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(linewidth = 0.25, color = "grey80"),
    panel.grid.minor = element_line(linewidth = 0.15, color = "grey92"),
    panel.spacing = unit(1, "lines"),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "plain", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

# ---- Kaydet ----------------------------------------------
out_dir <- file.path("outputs", "cross_scenario")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "fig_decade_progression_3ssp.png"), p,
       width = 13, height = 5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig_decade_progression_3ssp.pdf"), p,
       width = 13, height = 5, bg = "white")   # vektörel — Word'e gömmek için

message("✓ Grafik kaydedildi: ", file.path(out_dir, "fig_decade_progression_3ssp.png/.pdf"))
print(p)