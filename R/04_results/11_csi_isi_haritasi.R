## =========================================================
## 11) CSI Isı Haritası ve Trend  —  Ayrı Türkçe Script
## Grafikler tamamen TÜRKÇE, dosyalar _tr etiketli
## SSP1-2.6 (ssp126) ve SSP5-8.5 (ssp585) için
## Her senaryo KENDİ fig alt klasörüne kaydeder
## =========================================================
##
## GEREKSİNİMLER:
##   • Kütüphaneler         : dplyr, ggplot2, forcats, readr
##   • Brière fonksiyonları : briere_fn(), quad_lf_fn()
##   • Parametreler         : C_a, T0_a, Tm_a, C_eip, T0_e, Tm_e,
##                             C_lf, T0_lf, Tm_lf, lf_floor
##   • Tema & renkler       : theme_thesis(), COL_DISTRICT, DISTRICT_LABELS
##   • Yollar               : DIR_FIG (KÖK fig klasörü), DIR_TBL
##
## SENARYO VERİSİ (aşağıdaki YOL 1/2/3'ten biri yeterli — bkz. Sürücü):
##   Her senaryonun aylık iklim çerçevesi şu sütunları taşımalı:
##   temp_c, district_id, district_label, month_name, year
## =========================================================

library(dplyr)
library(ggplot2)
library(forcats)
library(readr)

## =========================================================
## Tek bir senaryo için CSI grafiklerini üreten fonksiyon
## =========================================================
csi_uret <- function(monthly, ssp_label, dir_fig, dir_tbl) {

  cat(">>> 11) CSI ısı haritası ve trend —", ssp_label, "<<<\n")

  ## --- Senaryoya özel fig klasörü (_tr çıktıları buraya) ---
  ssp_slug    <- gsub("[^A-Za-z0-9._-]", "_", ssp_label)   # klasör-güvenli ad
  DIR_FIG_SSP <- file.path(dir_fig, ssp_slug)
  dir.create(DIR_FIG_SSP, recursive = TRUE, showWarnings = FALSE)

  ## --- Brière fonksiyonunun teorik maksimumu ---
  ## Kapalı form karmaşık; geniş T aralığında sayısal maksimum alınır
  T_grid <- seq(0, 45, by = 0.1)
  a_max_theory   <- max(briere_fn(T_grid, C_a, T0_a, Tm_a))
  eip_max_theory <- max(briere_fn(T_grid, C_eip, T0_e, Tm_e))
  lf_max_theory  <- max(quad_lf_fn(T_grid, C_lf, T0_lf, Tm_lf))

  ## --- Aylık CSI (normalize bileşenlerin ortalaması) ---
  monthly_csi <- monthly %>%
    mutate(
      a_norm   = pmax(briere_fn(temp_c, C_a, T0_a, Tm_a), 0) / a_max_theory,
      lf_norm  = pmax(quad_lf_fn(temp_c, C_lf, T0_lf, Tm_lf), lf_floor) / lf_max_theory,
      eip_norm = pmax(briere_fn(temp_c, C_eip, T0_e, Tm_e), 0) / eip_max_theory,
      CSI = (a_norm + lf_norm + eip_norm) / 3
    )

  ## --- İlçe × Ay ortalaması → ısı haritası verisi ---
  csi_heat <- monthly_csi %>%
    group_by(district_label, month_name) %>%
    summarise(CSI_mean = mean(CSI, na.rm = TRUE), .groups = "drop") %>%
    ## Ay sırasını sabitle (Türkçe kısaltmalar korunur)
    mutate(month_name = factor(month_name,
                               levels = c("Oca","Şub","Mar","Nis","May","Haz",
                                          "Tem","Ağu","Eyl","Eki","Kas","Ara")))

  ## --- CSI ısı haritası grafiği (TÜRKÇE) ---
  fig_csi_heat <- ggplot(csi_heat,
                         aes(x = month_name,
                             y = fct_rev(factor(district_label)),
                             fill = CSI_mean)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", CSI_mean)), size = 2.8, colour = "grey20") +
    scale_fill_gradientn(
      colours = c("#EFF3FF","#BDD7E7","#6BAED6","#2171B5","#08306B"),
      name    = "İklim Uygunluk\nİndeksi (CSI)",
      limits  = c(0, 1),
      breaks  = c(0, 0.25, 0.50, 0.75, 1.00),
      labels  = c("0,00", "0,25", "0,50", "0,75", "1,00"),
      guide   = guide_colorbar(
        barwidth       = 18,
        barheight      = 1.0,
        title.position = "top",
        title.hjust    = 0.5,
        label.hjust    = 0.5
      )
    ) +
    labs(
      x       = "Ay",
      y       = "İlçe",
      title   = paste("İklim Uygunluk İndeksi \u2014", ssp_label),
      caption = "CSI = (a_norm + lf_norm + eip_norm) / 3; Bri\u00e8re termal performans e\u011frileri, Mordecai 2017"
    ) +
    theme_thesis() +
    theme(
      panel.grid   = element_blank(),
      axis.ticks   = element_blank(),
      plot.caption = element_text(size = 8, colour = "grey40", hjust = 0,
                                  margin = margin(t = 6))
    )

  ggsave("fig_csi_heat_tr.png", fig_csi_heat, path = DIR_FIG_SSP,
         width = 9, height = 5, dpi = 300)

  ## --- Yıllık CSI trendi (ilçe bazında doğrusal eğim) ---
  csi_yearly <- monthly_csi %>%
    group_by(district_id, district_label, year) %>%
    summarise(CSI_year = mean(CSI, na.rm = TRUE), .groups = "drop")

  fig_csi_trend <- ggplot(csi_yearly,
                          aes(x = year, y = CSI_year, colour = district_id)) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
    scale_colour_manual(
      values = COL_DISTRICT,
      labels = DISTRICT_LABELS,
      name   = NULL
    ) +
    scale_fill_manual(
      values = COL_DISTRICT,
      labels = DISTRICT_LABELS,
      guide  = "none"
    ) +
    labs(x = "Yıl", y = "Ortalama CSI",
         title = paste("CSI yıllık trend \u2014", ssp_label)) +
    theme_thesis()

  ggsave("fig_csi_trend_tr.png", fig_csi_trend, path = DIR_FIG_SSP,
         width = 8, height = 5, dpi = 300)

  ## --- Tablo çıktısı (senaryo etiketli, _tr) ---
  write_csv(csi_heat,
            file.path(dir_tbl, paste0("tbl_csi_monthly_", ssp_slug, "_tr.csv")))

  cat(">>> 11) tamamlandı →", DIR_FIG_SSP,
      "(fig_csi_heat_tr.png, fig_csi_trend_tr.png)\n")

  invisible(list(heat = csi_heat, yearly = csi_yearly))
}

## =========================================================
## SÜRÜCÜ  —  SSP1-2.6 ve SSP5-8.5
## ---------------------------------------------------------
## Etiket (grafikte görünen)  ->  kod (klasör/verideki karşılık)
SSP_HEDEF <- c("SSP1-2.6" = "ssp126",
               "SSP5-8.5" = "ssp585")

## Senaryonun aylık verisini getiren yardımcı.
## Ortamınızda VAR OLAN kaynağa göre otomatik seçer:
##   YOL 1) monthly_all : içinde 'ssp' sütunu olan tek birleşik çerçeve
##   YOL 2) monthly_ssp126 / monthly_ssp585 : senaryoya özel çerçeveler
get_monthly <- function(kod) {
  if (exists("monthly_all")) {
    df <- dplyr::filter(monthly_all, ssp == kod)   # 'ssp' sütun adınızsa
    if (nrow(df) > 0) return(df)
  }
  nesne <- paste0("monthly_", kod)                 # ör. monthly_ssp126
  if (exists(nesne)) return(get(nesne))
  stop(sprintf(
    "Senaryo verisi yok: '%s'. Ya 'ssp' sütunlu monthly_all ya da %s tanımlayın.",
    kod, nesne))
}

## Çoklu senaryo verisi mevcut mu?
coklu_var <- exists("monthly_all") ||
  (exists("monthly_ssp126") && exists("monthly_ssp585"))

if (coklu_var) {
  ## Her iki senaryoyu üret
  for (etiket in names(SSP_HEDEF)) {
    csi_uret(monthly   = get_monthly(SSP_HEDEF[[etiket]]),
             ssp_label = etiket,
             dir_fig   = DIR_FIG,
             dir_tbl   = DIR_TBL)
  }
} else if (exists("monthly") && exists("SSP_LABEL")) {
  ## GERİYE DÖNÜK UYUM: pipeline zaten aktif senaryoyu kurmuşsa
  ## (monthly + SSP_LABEL) yalnızca onu çalıştır — çökme yok.
  message("ℹ  Çoklu senaryo verisi bulunamadı; yalnızca aktif senaryo: ", SSP_LABEL)
  csi_uret(monthly, SSP_LABEL, DIR_FIG, DIR_TBL)
} else {
  stop(paste0(
    "Senaryo verisi bulunamadı. Şunlardan birini sağlayın:\n",
    "  • 'ssp' sütunlu monthly_all (ssp126/ssp585 değerleriyle), veya\n",
    "  • monthly_ssp126 ve monthly_ssp585 çerçeveleri, veya\n",
    "  • aktif senaryo için monthly + SSP_LABEL."))
}
