## =========================================================
## R/04_results/01_generate_ssp_outputs.R
## Her SSP senaryosu için grafik ve tablo üretimi
##
## Kullanım:
##   Sys.setenv(SSP_SCENARIO = "ssp245")
##   source("R/04_results/01_generate_ssp_outputs.R")
##
## Veya tüm SSP'ler için:
##   for (s in c("ssp126","ssp245","ssp585"))  {
##     Sys.setenv(SSP_SCENARIO = s)
##     source("R/04_results/01_generate_ssp_outputs.R")
##   }
##
## Çıktılar: outputs/{ssp}/figures/ ve outputs/{ssp}/tables/
## =========================================================



source("R/04_results/00_results_setup.R")

SSP <- Sys.getenv("SSP_SCENARIO", unset = "ssp245")
SSP_LABEL <- SSP_LABELS[SSP]

cat("\n", strrep("=", 55), "\n")
cat("Generating outputs for:", SSP_LABEL, "\n")
cat(strrep("=", 55), "\n\n")

# ---- Dizinler ----
DIR_FIG   <- here("outputs", SSP, "figures")
DIR_TBL   <- here("outputs", SSP, "tables")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TBL, recursive = TRUE, showWarnings = FALSE)

# ---- Veri yükle ----
dat <- load_ssp_data(SSP)
monthly <- dat$monthly
yearly  <- dat$yearly
horizon <- dat$horizon

gamma_val <- 1 / 5
TAU_VAL   <- 30L

cat("Data loaded:\n")
cat("  Monthly rows:", nrow(monthly), "\n")
cat("  Year range:", range(monthly$year), "\n")
cat("  Districts:", nlevels(monthly$district_id), "\n\n")

## =========================================================
## 1) İklim Koşulları
## =========================================================
p_temp_rh <- monthly %>%
  group_by(district_id, district_label, month_name) %>%
  summarise(T_mean  = mean(temp_c, na.rm = TRUE),
            RH_mean = mean(rh,     na.rm = TRUE),
            .groups = "drop") %>%
  ggplot(aes(x = month_name, group = district_id, colour = district_id)) +
  # 1. Sıcaklık Çizgileri (linetype aes içinde kalmalı)
  geom_line(aes(y = T_mean, linetype = "Sıcaklık (°C)"), linewidth = 0.8) +
  geom_point(aes(y = T_mean), size = 1.5) +
  # 2. Nem Çizgileri (Burada dışarıdaki linetype="dashed" kaldırıldı)
  geom_line(aes(y = RH_mean / 2.5, linetype = "Bağıl Nem (%)"), alpha = 0.6) + 
  scale_colour_manual(
    values = COL_DISTRICT,
    labels = DISTRICT_LABELS,
    name   = "İlçeler" # Karmaşayı önlemek için isim verilebilir
  ) +
  # Linetype manuel eşleştirme
  scale_linetype_manual(
    values = c("Sıcaklık (°C)" = "solid", "Bağıl Nem (%)" = "dashed"), 
    name = "Parametre"
  ) +
  scale_y_continuous(
    name = "Ortalama Sıcaklık (°C)",
    sec.axis = sec_axis(~ . * 2.5, name = "Bağıl Nem (%)")
  ) +
  labs(x     = NULL,
       title = paste("Mevsimsel Sıcaklık ve Nem Profili —", SSP_LABEL)) +
  theme_thesis() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",      # İlçeler ve Parametre üst üste binsin (yer tasarrufu)
    legend.margin = margin(t = 0)
  )

# Kaydetme (Word için ideal boyutlar)
ggsave("fig_temp_rh.png", p_temp_rh, path = DIR_FIG,
       width = 8.5, height = 5.5, dpi = 300)

## =========================================================
# 2) λ_local Isı Haritası
## =========================================================


cat(">>> 2) λ_local ısı haritası <<<\n")

heatmap_data <- monthly %>%
  group_by(district_label, month_name) %>%
  summarise(lam = mean(lambda_local_i1_mean, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(lam_plot = ifelse(lam < 1e-8, NA_real_, lam))

fig_heatmap <- ggplot(heatmap_data, aes(x = month_name, y = district_label,
                                        fill = lam_plot)) +
  geom_tile(colour = "white") +
  scale_fill_gradientn(
    colours  = c("grey92", "#FFF5B1", "#FEB24C", "#FC4E2A", "#B10026"),
    trans    = "log10",
    na.value = "grey88",
    breaks   = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
    labels   = label_scientific(),
    name     = expression(lambda[local]~"(log"[10]*" \u00f6l\u00e7ek)"),
    guide    = guide_colorbar(
      barwidth       = 20,
      barheight      = 1.2,
      title.position = "top",
      title.hjust    = 0.5
    )
  ) +
  labs(
    x       = NULL,
    y       = NULL,
    title   = paste("Ayl\u0131k yerel (otokton) bula\u015f h\u0131z\u0131 \u2014", SSP_LABEL),
    caption = paste0(
      "\u25a0 Koyu gri: Termal e\u015fik alt\u0131 (\u03bb < 10\u207b\u2078) \u2014 otokton bula\u015f imk\u00e2ns\u0131z\n",
      "\u25a0 A\u00e7\u0131k sar\u0131: Marjinal sezon ba\u015f\u0131/sonu (10\u207b\u2075\u201310\u207b\u00b3)\n",
      "\u25a0 Turuncu\u2013k\u0131rm\u0131z\u0131: Aktif bula\u015f sezonu (\u03bb \u2265 10\u207b\u00b2)"
    )
  ) +
  theme_thesis() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(
      size    = 8.5,
      colour  = "grey30",
      hjust   = 0,
      margin  = margin(t = 10)
    ),
    legend.title = element_text(size = 9, hjust = 0.5)
  )

ggsave("fig_heatmap_lambda.png", fig_heatmap, path = DIR_FIG,
       width = 9, height = 5.5, dpi = 300)


## =========================================================
## 3) P_est Mevsimselliği
## =========================================================
cat(">>> 3) P_est mevsimselliği <<<\n")

pest_season <- monthly %>%
  filter(p_establishment_mean > 0) %>%
  group_by(district_id, month_name) %>%
  summarise(pest_mean = mean(p_establishment_mean, na.rm = TRUE),
            .groups = "drop")

fig_pest_season <- ggplot(pest_season,
                          aes(x = month_name, y = pest_mean,
                              group = district_id,
                              colour = district_id)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_y_log10(labels = label_scientific()) +
  scale_colour_manual(values = COL_DISTRICT, name = NULL) +
  labs(x = NULL, y = expression(P[est]~"(log ölçek)"),
       title = paste("Yerleşme olasılığı mevsimsel profili —", SSP_LABEL)) +
  theme_thesis()

ggsave("fig_pest_season.png", fig_pest_season, path = DIR_FIG,
       width = 8, height = 5, dpi = 300)

## =========================================================
## 4) Yıllık Risk Trendi
## =========================================================
cat(">>> 4) Yıllık risk trendi <<<\n")

fig_yearly <- ggplot(yearly, aes(x = year, y = p_ge1_major_year_mean,
                                 colour = district_id,
                                 fill   = district_id)) +
  geom_ribbon(aes(ymin = p_ge1_major_year_p2_5,
                  ymax = p_ge1_major_year_p97_5),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(labels = label_scientific()) +
  scale_colour_manual(values = COL_DISTRICT, labels = DISTRICT_LABELS, name = NULL) +
  scale_fill_manual(values   = COL_DISTRICT, labels = DISTRICT_LABELS, guide = "none") +
  labs(x = "Yıl", y = expression(P["≥1 majör/yıl"]~"(log)"),
       title = paste("Yıllık salgın riski —", SSP_LABEL)) +
  theme_thesis()

ggsave("fig_yearly_risk.png", fig_yearly, path = DIR_FIG,
       width = 8, height = 5, dpi = 300)

## =========================================================
## 4) Annual Risk Trend  (ENGLISH)
## 01_generate_ssp_outputs.R içindeki "4) Yıllık Risk Trendi"
## bölümünün İngilizce karşılığı. Per-SSP döngü içinde
## (SSP_LABEL, DIR_FIG, yearly, COL_DISTRICT, DISTRICT_LABELS
## tanımlıyken) çalışır. Çıktı: fig_yearly_risk_en.png
##
## Yalnızca İngilizce isterseniz mevcut Türkçe bloğun yerine koyun;
## ikisini birlikte üretmek isterseniz bu bloğu Türkçe bloğun
## hemen altına ekleyin.
## =========================================================
cat(">>> 4) Annual risk trend <<<\n")

# İngilizce ilçe etiketleri: mevcut DISTRICT_LABELS'tan yalnızca İstanbul -> Istanbul
DISTRICT_LABELS_EN <- gsub("\u0130stanbul", "Istanbul", DISTRICT_LABELS)

fig_yearly_en <- ggplot(yearly, aes(x = year, y = p_ge1_major_year_mean,
                                    colour = district_id,
                                    fill   = district_id)) +
  geom_ribbon(aes(ymin = p_ge1_major_year_p2_5,
                  ymax = p_ge1_major_year_p97_5),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(labels = label_scientific()) +
  scale_colour_manual(values = COL_DISTRICT, labels = DISTRICT_LABELS_EN, name = NULL) +
  scale_fill_manual(values   = COL_DISTRICT, labels = DISTRICT_LABELS_EN, guide = "none") +
  labs(x = "Year",
       y = expression(P["\u22651 major/year"] ~ "(log)"),
       title = paste("Annual outbreak risk \u2014", SSP_LABEL)) +
  theme_thesis()

ggsave("fig_yearly_risk_en.png", fig_yearly_en, path = DIR_FIG,
       width = 8, height = 5, dpi = 300)

## =========================================================
## Not: Bu grafiği Rmd'de PNG olarak dahil ediyorsanız (ör. CSI trend
## chunk'ındaki gibi outputs/{SSP}/figures/... okuma), dosya yolunu
## "fig_yearly_risk_en.png" olarak güncelleyin.
## =========================================================

## =========================================================
## 5) Ufuk Tablosu
## =========================================================
cat(">>> 5) Ufuk risk tablosu <<<\n")

horizon_tbl <- horizon %>%
  select(district_label,
         p_mean = p_ge1_major_mean,
         p_lo = p_ge1_major_p2_5,
         p_hi = p_ge1_major_p97_5,
         Lambda = Lambda_import) %>%
  arrange(desc(p_mean)) %>%
  mutate(across(c(p_mean, p_lo, p_hi), ~formatC(.x, format = "e", digits = 2)),
         Lambda = round(Lambda, 1))

write_csv(horizon_tbl, file.path(DIR_TBL, "tbl_horizon.csv"))

ft <- flextable(horizon_tbl) %>% autofit() %>% theme_vanilla()
doc <- read_docx() %>% body_add_flextable(ft)
print(doc, target = file.path(DIR_TBL, "tbl_horizon.docx"))

# =========================================================
## 5b) İthalat Baskısı Özeti (λ_import trendi + kümülatif Λ)
## =========================================================
cat(">>> 5b) İthalat baskısı özeti <<<\n")

# ---- 5b.1) Yıllık ithalat baskısı trendi (log ölçek) ----
fig_import_trend <- ggplot(yearly,
                           aes(x = year, y = Lambda_import_year,
                               colour = district_id)) +
  geom_line(linewidth = 0.7) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.5,
              linetype = "dashed", alpha = 0.5) +
  scale_y_log10(labels = label_scientific()) +
  scale_colour_manual(values = COL_DISTRICT,
                      labels = DISTRICT_LABELS,
                      name = NULL) +
  labs(x = "Yıl",
       y = expression(Lambda["ithalat, yıllık"]~"(log ölçek)"),
       title = paste("Yıllık ithalat baskısı —", SSP_LABEL)) +
  theme_thesis()

ggsave("fig_import_trend.png", fig_import_trend, path = DIR_FIG,
       width = 8, height = 5, dpi = 300)

# ---- 5b.2) Aylık ithalat baskısı ısı haritası ----
import_heat <- monthly %>%
  group_by(district_label, month_name) %>%
  summarise(lambda_mean = mean(lambda_import, na.rm = TRUE),
            .groups = "drop")

fig_import_heat <- ggplot(import_heat,
                          aes(x = month_name, y = district_label,
                              fill = lambda_mean)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", trans = "log10",
                       labels = label_scientific(),
                       name = expression(lambda["ithalat"])) +
  labs(x = NULL, y = NULL,
       title = paste("Aylık ortalama ithalat baskısı —", SSP_LABEL)) +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("fig_import_heatmap.png", fig_import_heat, path = DIR_FIG,
       width = 8, height = 4, dpi = 300)

# ---- 5b.3) Kümülatif Λ_import özet tablosu ----
import_summary <- horizon %>%
  select(district_label, Lambda_import) %>%
  arrange(desc(Lambda_import)) %>%
  mutate(
    Lambda_fmt = round(Lambda_import, 1),
    Yorum = case_when(
      Lambda_import >= 100 ~ "Yüksek ithalat baskısı",
      Lambda_import >= 10  ~ "Orta ithalat baskısı",
      Lambda_import >= 1   ~ "Düşük ithalat baskısı",
      TRUE                 ~ "Çok düşük (< 1 beklenen vaka/50 yıl)"
    )
  ) %>%
  select(`İlçe` = district_label,
         `Λ_ithalat (50 yıl)` = Lambda_fmt,
         Yorum)

write_csv(import_summary, file.path(DIR_TBL, "tbl_import_summary.csv"))

# Flextable ile docx
ft_imp <- flextable(import_summary) %>% autofit() %>% theme_vanilla()
doc_imp <- read_docx() %>% body_add_flextable(ft_imp)
print(doc_imp, target = file.path(DIR_TBL, "tbl_import_summary.docx"))

cat("  İthalat baskısı grafikleri ve tablosu kaydedildi.\n")


## =========================================================
## 6) Model Doğrulaması
## =========================================================
cat(">>> 6) Model doğrulaması <<<\n")

val_df <- monthly %>%
  mutate(
    rho_ratio = ifelse(lambda_local_i1_mean > 0,
                       gamma_val / lambda_local_i1_mean, Inf),
    P_est_analytic = case_when(
      lambda_local_i1_mean <= 0  ~ 0,
      abs(rho_ratio - 1) < 1e-10 ~ 1 / TAU_VAL,
      TRUE ~ (1 - rho_ratio) / (1 - rho_ratio^TAU_VAL)
    ),
    abs_diff = abs(p_establishment_mean - P_est_analytic),
    active   = lambda_local_i1_mean > 0
  )

val_active <- filter(val_df, active)

val_summary <- tibble(
  Metrik = c("Toplam kombinasyon", "Aktif (λ>0)",
             "Ort. |fark|", "Maks |fark|", "< 0.02 eşiği"),
  Değer = c(nrow(val_df), nrow(val_active),
            formatC(mean(val_active$abs_diff), format = "e", digits = 2),
            formatC(max(val_active$abs_diff), format = "e", digits = 2),
            paste0(sum(val_active$abs_diff < 0.02), "/", nrow(val_active)))
)

write_csv(val_summary, file.path(DIR_TBL, "tbl_validation.csv"))

fig_val <- ggplot(val_active, aes(x = P_est_analytic,
                                   y = p_establishment_mean,
                                   colour = district_label)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  scale_x_log10(labels = label_scientific()) +
  scale_y_log10(labels = label_scientific()) +
  scale_colour_manual(values = COL_DISTRICT, name = NULL) +
  labs(x = expression(P[est]^{analitik}),
       y = expression(P[est]^{MC}),
       title = paste("Model doğrulama —", SSP_LABEL)) +
  theme_thesis()

ggsave("fig_validation.png", fig_val, path = DIR_FIG,
       width = 7, height = 5, dpi = 300)

## =========================================================
## 7) İthalat Duyarlılık (k + η)
## =========================================================
cat(">>> 7) İthalat duyarlılığı <<<\n")

sens_imp <- load_sens_imp(SSP)
if (!is.null(sens_imp)) {
  write_csv(sens_imp, file.path(DIR_TBL, "tbl_sens_importation.csv"))
  cat("  İthalat duyarlılık tablosu kaydedildi.\n")
} else {
  cat("  İthalat duyarlılık verisi bulunamadı — atlanıyor.\n")
}

## =========================================================
## 8) Model Duyarlılık (m, beta_vh, IP)
## =========================================================
cat(">>> 8) Model duyarlılığı <<<\n")

sens_mc <- load_sens_mc(SSP)
if (!is.null(sens_mc)) {
  # Tornado plot — signed log10 scale
  sens_tor <- read_csv(file.path(here("outputs", SSP, "sensitivity", "ctmc_mc"),
                                 "sensitivity_tornado.csv"),
                       show_col_types = FALSE) %>%
    mutate(
      # Signed log10 transform: sign(x) * log10(|x| + 1)
      signed_log10 = sign(mean_delta_pct) * log10(abs(mean_delta_pct) + 1),
      # Turkish scenario labels
      senaryo_tr = dplyr::recode(scenario, !!!SCEN_TR)
    )

  # Symmetric axis limits
  max_abs <- max(abs(sens_tor$signed_log10), na.rm = TRUE) * 1.1

  fig_tornado <- ggplot(sens_tor,
                        aes(x = reorder(senaryo_tr, abs(signed_log10)),
                            y = signed_log10)) +
    geom_col(aes(fill = signed_log10 > 0), show.legend = FALSE, width = 0.7) +
    geom_text(aes(label = sprintf("%+.0f%%", mean_delta_pct),
                  hjust = ifelse(signed_log10 >= 0, -0.1, 1.1)),
              size = 3, colour = "grey30") +
    coord_flip(ylim = c(-max_abs, max_abs)) +
    scale_fill_manual(values = c("TRUE" = "#E63946", "FALSE" = "#457B9D")) +
    labs(x = NULL,
         y = expression("Baz senaryoya göre değişim — sign(x) " %*% " log"[10]*"(|%Δ|+1)"),
         title = paste("Parametre duyarlılığı —", SSP_LABEL),
         caption = "Çubuk yanındaki etiketler orijinal % değişimi gösterir") +
    theme_thesis() +
    theme(axis.text.y = element_text(size = 9))

  ggsave("fig_tornado.png", fig_tornado, path = DIR_FIG,
         width = 8, height = 5, dpi = 300)
  cat("  Tornado grafiği (log10) kaydedildi.\n")
} else {
  cat("  Model duyarlılık verisi bulunamadı — atlanıyor.\n")
}

## =========================================================
## 9) GBD Ülke Katkıları
## =========================================================
cat(">>> 9) GBD ülke katkıları <<<\n")

cc <- load_country_contrib(SSP)
if (!is.null(cc)) {
  top_25 <- cc %>% slice_max(contribution_pct, n =25)

  fig_country <- ggplot(top_25,
                        aes(x = reorder(country, contribution_pct),
                            y = contribution_pct)) +
    geom_col(fill = "#E63946", alpha = 0.85) +
    coord_flip() +
    labs(x = NULL, y = "İthalat riskine katkı (%)",
         title = "Kaynak ülke bazlı dang ithalat riski",
         subtitle = paste("GBD 2023 × turist ağırlığı —", SSP_LABEL)) +
    theme_thesis()

  ggsave("fig_country_contrib.png", fig_country, path = DIR_FIG,
         width = 7, height = 5, dpi = 300)

  write_csv(top_25, file.path(DIR_TBL, "tbl_country_contrib.csv"))
  cat("  Ülke katkısı grafiği ve tablosu kaydedildi.\n")
} else {
  cat("  Ülke katkısı verisi bulunamadı — atlanıyor.\n")
}

## =========================================================
## 10) LHS–PRCC Duyarlılık Analizi (SSP-spesifik iklim verisi)
## =========================================================
cat(">>> 10) LHS–PRCC duyarlılık analizi (SSP-spesifik) <<<\n")

library(lhs)

set.seed(123)
n_lhs <- 2000

# ---- SSP-spesifik: T ve RH değerlerini gerçek iklim verisinden örnekle ----
# Aktif bulaş aylarından (lambda_local > 0) rasgele örnekleme
active_climate <- monthly %>%
  filter(lambda_local_i1_mean > 0) %>%
  select(temp_c, rh) %>%
  filter(is.finite(temp_c), is.finite(rh))

if (nrow(active_climate) < 100) {
  active_climate <- monthly %>%
    select(temp_c, rh) %>%
    filter(is.finite(temp_c), is.finite(rh))
}

# Gerçek iklim verisinden bootstrap örnekleme (sıralama korunmaz → PRCC farklılaşır)
climate_sample_idx <- sample(nrow(active_climate), n_lhs, replace = TRUE)
T_sampled  <- active_climate$temp_c[climate_sample_idx]
RH_sampled <- active_climate$rh[climate_sample_idx]

T_range  <- round(range(T_sampled), 1)
RH_range <- round(range(RH_sampled), 1)

cat("  SSP-spesifik iklim (bootstrap örnekleme):\n")
cat("    T_C  :", T_range[1], "–", T_range[2], "°C\n")
cat("    RH   :", RH_range[1], "–", RH_range[2], "%\n")
cat("    Aktif ay havuzu:", nrow(active_climate), "gözlem\n")

# ---- Diğer parametreler: LHS ile uniform örnekleme ----
param_ranges_other <- list(
  m       = c(0.10, 2.00),
  beta_vh = c(0.10, 0.60),
  beta_hv = c(0.10, 0.60),
  ip_days = c(3, 10)
)

lhs_unit <- randomLHS(n_lhs, length(param_ranges_other))
colnames(lhs_unit) <- names(param_ranges_other)

lhs_params <- as.data.frame(lhs_unit)
for (nm in names(param_ranges_other)) {
  rng <- param_ranges_other[[nm]]
  lhs_params[[nm]] <- rng[1] + lhs_unit[, nm] * (rng[2] - rng[1])
}

# T_C ve RH: gerçek iklim verisinden örneklenmiş değerler
lhs_params$T_C <- T_sampled
lhs_params$RH  <- RH_sampled

# Mordecai 2017 Ae. albopictus parametreleri
C_a <- 1.93e-4; T0_a <- 10.25; Tm_a <- 38.32
C_eip <- 1.09e-4; T0_e <- 10.39; Tm_e <- 43.05
C_lf <- 1.43e-1; T0_lf <- 6.24; Tm_lf <- 38.25
k_vpd <- 0.5; SVPD_ref <- 1.0; lf_floor <- 0.25

briere_fn <- function(T, c, T0, Tm) ifelse(T > T0 & T < Tm, c * T * (T - T0) * sqrt(Tm - T), 0)
quad_lf_fn <- function(T, c, T0, Tm) pmax(c * (T - T0) * (Tm - T), lf_floor)
svpd_fn <- function(T, RH) pmax(0.6108 * exp(17.27 * T / (T + 237.3)) * (1 - RH / 100), 0)

compute_R0_lhs <- function(p) {
  a <- briere_fn(p$T_C, C_a, T0_a, Tm_a)
  lf <- quad_lf_fn(p$T_C, C_lf, T0_lf, Tm_lf)
  eip <- 1 / pmax(briere_fn(p$T_C, C_eip, T0_e, Tm_e), 1e-6)
  vpd <- svpd_fn(p$T_C, p$RH)
  mu_v <- (1 / lf) * exp(k_vpd * (vpd - SVPD_ref))
  gamma <- 1 / p$ip_days
  pmax(p$m * a^2 * p$beta_vh * p$beta_hv * exp(-mu_v * eip) / (mu_v * gamma), 0)
}

lhs_params$R0 <- apply(lhs_params, 1, function(row) compute_R0_lhs(as.list(row)))

param_names <- c("m", "beta_vh", "beta_hv", "T_C", "RH", "ip_days")
prcc_result <- purrr::map_dfr(param_names, function(nm) {
  others <- setdiff(param_names, nm)
  rank_xi <- rank(lhs_params[[nm]])
  rank_y <- rank(lhs_params$R0)
  rank_others <- lhs_params[, others] %>% mutate(across(everything(), rank))
  resid_xi <- residuals(lm(rank_xi ~ ., data = rank_others))
  resid_y <- residuals(lm(rank_y ~ ., data = rank_others))
  ct <- cor.test(resid_xi, resid_y, method = "pearson")
  tibble(parameter = nm, PRCC = ct$estimate, p_value = ct$p.value,
         ci_lo = ct$conf.int[1], ci_hi = ct$conf.int[2])
})

write_csv(prcc_result, file.path(DIR_TBL, "tbl_prcc.csv"))

# Parametre aralıklarını kaydet
param_ranges_tbl <- tibble(
  Parameter = c("m", "beta_vh", "beta_hv", "ip_days", "T_C", "RH"),
  Lower = c(
    param_ranges_other$m[1], 
    param_ranges_other$beta_vh[1], 
    param_ranges_other$beta_hv[1], 
    param_ranges_other$ip_days[1], 
    T_range[1], 
    RH_range[1]
  ),
  Upper = c(
    param_ranges_other$m[2], 
    param_ranges_other$beta_vh[2], 
    param_ranges_other$beta_hv[2], 
    param_ranges_other$ip_days[2], 
    T_range[2], 
    RH_range[2]
  ),
  Source = c(
    "Literature", "Literature", "Literature", "Literature",
    paste0("Bootstrap (", SSP_LABEL, " active months)"),
    paste0("Bootstrap (", SSP_LABEL, " active months)")
  )
)

# Tabloyu kaydet
write_csv(param_ranges_tbl, file.path(DIR_TBL, "tbl_param_ranges.csv"))

fig_prcc <- ggplot(prcc_result, aes(x = reorder(parameter, abs(PRCC)), y = PRCC)) +
  geom_col(aes(fill = PRCC > 0), show.legend = FALSE, width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E63946", "FALSE" = "#457B9D")) +
  labs(x = NULL, y = "PRCC",
       title = paste(SSP_LABEL, "LHS\u2013PRCC Sensitivity Analysis"),
       subtitle = sprintf("n = %d | T: %.1f\u2013%.1f\u00b0C | RH: %.1f\u2013%.1f%% (active-month bootstrap)",
                          n_lhs, T_range[1], T_range[2],
                          RH_range[1], RH_range[2])) +
  theme_thesis()

ggsave("fig_prcc.png", fig_prcc, path = DIR_FIG, width = 7, height = 4.5, dpi = 300)

## =========================================================
## 11) CSI Isı Haritası ve Trend
## =========================================================
cat(">>> 11) CSI ısı haritası ve trend <<<\n")

# Brière fonksiyonunun teorik maksimumu
# T_opt = (2*Tm*T0 + Tm^2 - T0^2 ± ...) / 3  — kapalı form karmaşık
# Pratik: geniş T aralığında sayısal maksimum
T_grid <- seq(0, 45, by = 0.1)
a_max_theory   <- max(briere_fn(T_grid, C_a, T0_a, Tm_a))
eip_max_theory <- max(briere_fn(T_grid, C_eip, T0_e, Tm_e))
lf_max_theory  <- max(quad_lf_fn(T_grid, C_lf, T0_lf, Tm_lf))

monthly_csi <- monthly %>%
  mutate(
    a_norm   = pmax(briere_fn(temp_c, C_a, T0_a, Tm_a), 0) / a_max_theory,
    lf_norm  = pmax(quad_lf_fn(temp_c, C_lf, T0_lf, Tm_lf), lf_floor) / lf_max_theory,
    eip_norm = pmax(briere_fn(temp_c, C_eip, T0_e, Tm_e), 0) / eip_max_theory,
    CSI = (a_norm + lf_norm + eip_norm) / 3
  )

csi_heat <- monthly_csi %>%
  group_by(district_label, month_name) %>%
  summarise(CSI_mean = mean(CSI, na.rm = TRUE), .groups = "drop")

csi_heat <- csi_heat %>%
  mutate(month_name = factor(month_name,
                             levels = c("Oca","Şub","Mar","Nis","May","Haz",
                                        "Tem","Ağu","Eyl","Eki","Kas","Ara"),
                             labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                        "Jul","Aug","Sep","Oct","Nov","Dec")))

fig_csi_heat <- ggplot(csi_heat, aes(x = month_name, y = fct_rev(factor(district_label)),
                                     fill = CSI_mean)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", CSI_mean)), size = 2.8, colour = "grey20") +
  scale_fill_gradientn(
    colours = c("#EFF3FF","#BDD7E7","#6BAED6","#2171B5","#08306B"),
    name    = "Climate Suitability\nIndex (CSI)",
    limits  = c(0, 1),
    breaks  = c(0, 0.25, 0.50, 0.75, 1.00),
    labels  = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    guide   = guide_colorbar(
      barwidth       = 18,
      barheight      = 1.0,
      title.position = "top",
      title.hjust    = 0.5,
      label.hjust    = 0.5
    )
  ) +
  labs(
    x       = "Year",
    y       = "Mean CSI",
    title   = paste("Climate Suitability Index \u2014", SSP_LABEL),
    caption = "CSI = (a_norm + lf_norm + eip_norm) / 3; Bri\u00e8re thermal performance curves, Mordecai 2017"
  ) +
  theme_thesis() +
  theme(
    panel.grid   = element_blank(),
    axis.ticks   = element_blank(),
    plot.caption = element_text(size = 8, colour = "grey40", hjust = 0,
                                margin = margin(t = 6))
  )

ggsave("fig_csi_heat.png", fig_csi_heat, path = DIR_FIG,
       width = 9, height = 5, dpi = 300)

csi_yearly <- monthly_csi %>%
  group_by(district_id, district_label, year) %>%
  summarise(CSI_year = mean(CSI, na.rm = TRUE), .groups = "drop")

fig_csi_trend <- ggplot(csi_yearly, aes(x = year, y = CSI_year,
                                        colour = district_id)) +  # ← TEK DEĞİŞİKLİK
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
  scale_colour_manual(
    values = COL_DISTRICT,
    labels = DISTRICT_LABELS,
    name = NULL
  ) +
  scale_fill_manual(
    values = COL_DISTRICT,
    labels = DISTRICT_LABELS,
    guide = "none"
  ) +
  labs(x = "Yıl", y = "Ortalama CSI",
       title = paste("CSI anuual trend —", SSP_LABEL)) +
  theme_thesis()

ggsave("fig_csi_trend.png", fig_csi_trend, path = DIR_FIG, width = 8, height = 5, dpi = 300)

write_csv(csi_heat, file.path(DIR_TBL, "tbl_csi_monthly.csv"))

## =========================================================
## 12) Moran's I Mekânsal Otokorelasyon
## =========================================================
cat(">>> 12) Moran's I <<<\n")

if (requireNamespace("spdep", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
  library(spdep); library(sf)

  district_coords <- data.frame(
    district_id = c("TUR.10.4_1","TUR.39.3_1","TUR.40.25_1","TUR.59.4_1","TUR.81.6_1"),
    label = c("Hopa","Eğirdir","Kartal","Fethiye","Zonguldak"),
    lon = c(41.12, 30.85, 29.19, 29.12, 31.79),
    lat = c(41.41, 37.88, 40.89, 36.62, 41.45)
  )

  pts_sf <- st_as_sf(district_coords, coords = c("lon","lat"), crs = 4326) %>%
    st_transform(crs = 32636)

  coords_mat <- st_coordinates(pts_sf)
  knn2 <- knearneigh(coords_mat, k = 2)
  nb2 <- knn2nb(knn2)
  w2 <- nb2listw(nb2, style = "W")

  horizon_vec <- horizon %>%
    arrange(factor(district_id, levels = district_coords$district_id)) %>%
    pull(p_ge1_major_mean)

  csi_vec <- monthly_csi %>%
    group_by(district_id) %>%
    summarise(CSI_mean = mean(CSI, na.rm = TRUE), .groups = "drop") %>%
    arrange(factor(district_id, levels = district_coords$district_id)) %>%
    pull(CSI_mean)

  mt_p <- moran.test(log10(pmax(horizon_vec, 1e-20)), w2,
                     randomisation = TRUE, alternative = "two.sided")
  mt_csi <- moran.test(csi_vec, w2, randomisation = TRUE, alternative = "two.sided")

  moran_tbl <- tibble(
    Değişken = c("log₁₀(p_ge1_major)", "CSI"),
    `Moran's I` = c(round(mt_p$estimate[1], 4), round(mt_csi$estimate[1], 4)),
    `p-değeri` = c(format.pval(mt_p$p.value, digits = 3),
                   format.pval(mt_csi$p.value, digits = 3)),
    Yorum = c(ifelse(mt_p$p.value < 0.05, "Anlamlı kümelenme", "Rastgele"),
              ifelse(mt_csi$p.value < 0.05, "Anlamlı kümelenme", "Rastgele"))
  )

  write_csv(moran_tbl, file.path(DIR_TBL, "tbl_moran.csv"))

  z_p <- scale(log10(pmax(horizon_vec, 1e-20)))[,1]
  z_csi <- scale(csi_vec)[,1]
  lag_p <- lag.listw(w2, z_p)
  lag_csi <- lag.listw(w2, z_csi)

  moran_df <- tibble(label = district_coords$label,
                     z_p = z_p, lag_p = lag_p,
                     z_csi = z_csi, lag_csi = lag_csi)

  p_moran <- ggplot(moran_df, aes(x = z_p, y = lag_p)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_smooth(method = "lm", se = TRUE, colour = "#E63946", alpha = 0.15) +
    geom_point(size = 4, colour = "#E63946", shape = 21, fill = "white", stroke = 1.5) +
    ggrepel::geom_text_repel(aes(label = label), size = 3.2) +
    labs(x = "Std. log₁₀(p_ge1_major)", y = "Mekânsal gecikme",
         title = paste("Moran saçılımı —", SSP_LABEL)) +
    theme_thesis()

  ggsave("fig_moran.png", p_moran, path = DIR_FIG, width = 6, height = 5, dpi = 300)
  cat("  Moran's I tamamlandı.\n")
} else {
  cat("  spdep/sf paketi yok — Moran's I atlanıyor.\n")
}

## =========================================================
## 13) Dekadal risk progresyonu
## =========================================================
cat(">>> 13) Dekadal risk progresyonu <<<\n")

decade_df <- yearly %>%
  mutate(
    dekad_yil = floor(year / 10) * 10,
    dekad = factor(
      dekad_yil,
      levels = seq(2020, 2070, 10),
      labels = c("2020'ler","2030'lar","2040'lar",
                 "2050'ler","2060'lar","2070'ler")
    )
  ) %>%
  filter(!is.na(dekad)) %>%
  group_by(district_id, district_label, dekad) %>%
  summarise(p_ort = mean(p_ge1_major_year_mean, na.rm = TRUE), .groups = "drop")

# Grafik
fig_decade <- ggplot(decade_df, aes(x = dekad, y = p_ort,
                                    colour = district_id,
                                    group  = district_id)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = COL_DISTRICT, labels = DISTRICT_LABELS, name = "\u0130l\u00e7e") +
  scale_y_log10(
    labels = label_scientific(digits = 1),
    breaks = 10^seq(-9, 0, 1)
  ) +
  labs(
    title = paste("On y\u0131ll\u0131k risk progresyonu \u2014", SSP_LABEL),
    x     = "On Y\u0131l",
    y     = "Ort. y\u0131ll\u0131k risk (log\u2081\u2080)"
  ) +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("fig_decade.png", fig_decade, path = DIR_FIG,
       width = 9, height = 5, dpi = 300)

# Tablo — ham sayısal değerlerle yaz
decade_wide <- decade_df %>%
  select(district_label, dekad, p_ort) %>%
  pivot_wider(names_from = dekad, values_from = p_ort) %>%
  rename("\u0130l\u00e7e" = district_label)

write_csv(decade_wide, file.path(DIR_TBL, "tbl_decade.csv"))
cat("  tbl_decade.csv yazıldı:", unique(yearly$ssp), "\n")



## =========================================================
## 14) Bulaş Sezonu Uzunluğu Değişimi
## =========================================================
cat(">>> 14) Bulaş sezonu uzunluğu <<<\n")

LAMBDA_THRESH <- 1e-4

season_yr <- monthly %>%
  group_by(district_id, district_label, year) %>%
  summarise(
    season_len = sum(lambda_local_i1_mean > LAMBDA_THRESH, na.rm = TRUE),
    peak_month = AY_TR[which.max(lambda_local_i1_mean)],
    .groups = "drop"
  )

fig_season <- ggplot(season_yr, aes(x = year, y = season_len,
                                     colour = district_id, fill = district_id)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.2) +
  scale_colour_manual(values = COL_DISTRICT, labels = DISTRICT_LABELS, name = "İlçe") +
  scale_fill_manual(values = COL_DISTRICT, labels = DISTRICT_LABELS, name = "İlçe") +
  scale_y_continuous(breaks = 0:8) +
  labs(title = paste("Bulaş sezonu uzunluğu —", SSP_LABEL),
       x = "Yıl", y = "Aktif ay sayısı (λ_local > 0)") +
  theme_thesis()

ggsave("fig_season.png", fig_season, path = DIR_FIG, width = 8, height = 5, dpi = 300)

season_comp <- season_yr %>%
  mutate(donem = case_when(
    year >= 2025 & year <= 2035 ~ "Erken",
    year >= 2065 & year <= 2075 ~ "Geç",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(donem)) %>%
  group_by(district_label, donem) %>%
  summarise(sezon_ort = mean(season_len, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = donem, values_from = sezon_ort) %>%
  mutate(
    Erken = round(Erken, 1),
    Geç   = round(Geç,   1),
    delta = round(Geç - Erken, 1)
  )

write_csv(season_comp, file.path(DIR_TBL, "tbl_season.csv"))



## =========================================================
## 15) Çoklu Regresyon (Aktif bulaş sezonu bazlı)
## =========================================================
cat(">>> 15) Çoklu regresyon (bulaş sezonu) <<<\n")

# ---- Yalnızca aktif aylardan yıllık istatistik türet ----
stat_df <- monthly %>%
  filter(lambda_local_i1_mean > 0) %>%
  group_by(district_id, district_label, year) %>%
  summarise(
    T_season    = mean(temp_c, na.rm = TRUE),
    RH_season   = mean(rh, na.rm = TRUE),
    n_active    = n(),
    pest_season = mean(p_establishment_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    yearly %>% select(district_id, year,
                      p_ge1_major_year_mean, Lambda_import_year),
    by = c("district_id", "year")
  ) %>%
  filter(p_ge1_major_year_mean > 0) %>%
  mutate(log_p = log10(p_ge1_major_year_mean))

if (nrow(stat_df) > 10) {
  lm1 <- lm(log_p ~ T_season + I(T_season^2), data = stat_df)
  lm2 <- lm(log_p ~ T_season + I(T_season^2) + RH_season, data = stat_df)
  lm3 <- lm(log_p ~ T_season + I(T_season^2) + RH_season +
              log10(Lambda_import_year + 1), data = stat_df)
  
  # Optimum sıcaklık hesabı
  b1 <- coef(lm3)["T_season"]
  b2 <- coef(lm3)["I(T_season^2)"]
  T_opt <- -b1 / (2 * b2)
  cat("  Bulaş sezonu optimum T:", round(T_opt, 1), "°C\n")
  
  reg_comp <- tibble(
    Model = c("M1: T + T²", "M2: + RH", "M3: + Λ_import"),
    R2 = c(summary(lm1)$r.squared, summary(lm2)$r.squared,
           summary(lm3)$r.squared),
    adj_R2 = c(summary(lm1)$adj.r.squared, summary(lm2)$adj.r.squared,
               summary(lm3)$adj.r.squared),
    delta_R2 = c(NA, summary(lm2)$r.squared - summary(lm1)$r.squared,
                 summary(lm3)$r.squared - summary(lm2)$r.squared),
    ANOVA_p = c(NA, anova(lm1, lm2)$`Pr(>F)`[2],
                anova(lm2, lm3)$`Pr(>F)`[2])
  ) %>%
    mutate(across(c(R2, adj_R2, delta_R2), ~round(.x, 4)),
           ANOVA_p = ifelse(is.na(ANOVA_p), "—",
                            format.pval(ANOVA_p, digits = 3)))
  
  write_csv(reg_comp, file.path(DIR_TBL, "tbl_regression.csv"))
  
  m3_coef <- broom::tidy(lm3, conf.int = TRUE) %>%
    transmute(
      Terim = case_when(
        term == "(Intercept)" ~ "Sabit",
        term == "T_season" ~ "T_sezon",
        term == "I(T_season^2)" ~ "T_sezon²",
        term == "RH_season" ~ "RH_sezon",
        term == "log10(Lambda_import_year + 1)" ~ "log₁₀(Λ+1)"
      ),
      Katsayı = round(estimate, 4),
      SE = round(std.error, 4),
      p = format.pval(p.value, digits = 3)
    )
  
  write_csv(m3_coef, file.path(DIR_TBL, "tbl_regression_coef.csv"))
  cat("  R² final model:", round(summary(lm3)$r.squared, 3), "\n")
  cat("  T_opt (sezon):", round(T_opt, 1), "°C\n")
} else {
  cat("  Yetersiz veri — regresyon atlanıyor.\n")
}

cat("\n", strrep("=", 55), "\n")
cat("DONE:", SSP_LABEL, "\n")
cat("Figures:", DIR_FIG, "\n")
cat("Tables:", DIR_TBL, "\n")
cat(strrep("=", 55), "\n")



## =========================================================
## 16) Linear Mixed Regression (Random effect= District)
## =========================================================

library(lme4)
library(broom.mixed)

stat_monthly <- monthly %>%
  filter(lambda_local_i1_mean > 0,
         p_month_major_mean > 0) %>%
  mutate(
    log_p   = log10(p_month_major_mean),
    log_lam = log10(lambda_import + 1e-10)
  )

# ---- Model 1: Random intercept (ilçe) + iklim ----
mm1 <- lmer(log_p ~ temp_c + I(temp_c^2) + rh +
              (1 | district_id),
            data = stat_monthly)

# ---- Model 2: + ithalat baskısı ----
mm2 <- lmer(log_p ~ temp_c + I(temp_c^2) + rh + log_lam +
              (1 | district_id),
            data = stat_monthly)

# ---- Model 3: Random slope (sıcaklık etkisi ilçeye göre değişir mi?) ----
mm3 <- lmer(log_p ~ temp_c + I(temp_c^2) + rh + log_lam +
              (1 + temp_c | district_id),
            data = stat_monthly)

# Optimum sıcaklık
b1 <- fixef(mm1)["temp_c"]
b2 <- fixef(mm1)["I(temp_c^2)"]
T_opt <- -b1 / (2 * b2)
cat("Mixed model T_opt:", round(T_opt, 1), "°C\n")

# Model karşılaştırma
anova(mm1, mm2, mm3)

# Katsayı tablosu
tidy(mm2, conf.int = TRUE, effects = "fixed")

# Random effects: ilçeler arası varyans
VarCorr(mm2)

# ICC: varyansın ne kadarı ilçeler arası
icc <- as.data.frame(VarCorr(mm2))$vcov[1] /
  (as.data.frame(VarCorr(mm2))$vcov[1] + sigma(mm2)^2)
cat("ICC:", round(icc, 3), "\n")


