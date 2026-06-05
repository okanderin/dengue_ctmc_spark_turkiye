## bulgular.Rmd — Revizyon Haritası
## Tarih: 2026-03-30
## Amaç: GBD-weighted importation, S_m relocasyonu, ×12 düzeltmesi,
##        SSP-aware çıktılar ve ρ kaldırılması sonrası gerekli değişiklikler

---

## DEĞİŞİKLİK 1 — Dosya Yolları (SSP-aware)  [ZORUNLU]

### Satır 119-121: Sabit yollar → SSP-aware yollar

ESKİ:
```r
DIR_SIM  <- here("outputs", "simulation")
DIR_SENS <- here("outputs", "sensitivity", "ctmc_mc")
DIR_PLOTS <- here("outputs", "figures")
DIR_TABLO <- here("outputs", "tablolar")
```

YENİ:
```r
SSP_SCENARIO <- Sys.getenv("SSP_SCENARIO", unset = "ssp245")
SSP_LABEL    <- toupper(gsub("ssp", "SSP", SSP_SCENARIO))

DIR_SSP_OUT  <- here("outputs", SSP_SCENARIO)
DIR_SIM      <- file.path(DIR_SSP_OUT, "simulation")
DIR_SENS_IMP <- file.path(DIR_SSP_OUT, "sensitivity", "importation")
DIR_SENS_MC  <- file.path(DIR_SSP_OUT, "sensitivity", "ctmc_mc")
DIR_PLOTS    <- file.path(DIR_SSP_OUT, "figures")
DIR_TABLO    <- file.path(DIR_SSP_OUT, "tables")
DIR_DIAG     <- file.path(DIR_SSP_OUT, "diagnostics")

dir.create(DIR_PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLO, recursive = TRUE, showWarnings = FALSE)
```

### Başlık (satır 2-3): SSP bilgisi ekle

YENİ:
```yaml
title: "Türkiye'de Dengue Yerleşme Riski: CTMC Spark Model Sonuçları"
subtitle: "`r paste0('Bulgular — ', toupper(gsub('ssp','SSP',Sys.getenv('SSP_SCENARIO','ssp245'))), ' Senaryosu (2025–2075)')`"
```

---

## DEĞİŞİKLİK 2 — İthalat Duyarlılık Verisi [ZORUNLU]

### Satır 146-170: Eski sensitivity yapısı → yeni (η-only + k)

Eski yapı: `sensitivity_summary.csv` (9 senaryo: m, beta_vh, ip_days)
Yeni yapı: İKİ ayrı duyarlılık kaynağı:
  (a) İthalat duyarlılığı: `DIR_SENS_IMP/sensitivity_summary_{ssp}.csv`
      → k × 3 düzey + η × 3 düzey
  (b) Model duyarlılığı: `DIR_SENS_MC/sensitivity_summary.csv`
      → m × 5 düzey + beta_vh ± 20% + ip_days ± 20%

YENİ KOD (satır 146 yerine):
```r
# ── İthalat duyarlılık verisi ──────────────────────────────────────
SENS_IMP_PATH <- file.path(DIR_SENS_IMP,
                           paste0("sensitivity_summary_", SSP_SCENARIO, ".csv"))
SENS_IMP_READY <- file.exists(SENS_IMP_PATH)

if (SENS_IMP_READY) {
  sens_imp <- read_csv(SENS_IMP_PATH, show_col_types = FALSE)
  cat("İthalat duyarlılık verisi yüklendi:", nrow(sens_imp), "satır\n")
} else {
  message("⚠ İthalat duyarlılık verisi bulunamadı — ilgili bölümler atlanacak.")
}

# ── Model (CTMC) duyarlılık verisi ─────────────────────────────────
SENS_MC_PATH <- file.path(DIR_SENS_MC, "sensitivity_summary.csv")
SENS_MC_READY <- file.exists(SENS_MC_PATH)

if (SENS_MC_READY) {
  sens_sum <- read_csv(SENS_MC_PATH, show_col_types = FALSE) %>%
    mutate(
      district_id    = factor(district_id, levels = DISTRICT_ORDER),
      district_label = DISTRICT_LABELS[as.character(district_id)]
    )
  sens_tor <- read_csv(file.path(DIR_SENS_MC, "sensitivity_tornado.csv"),
                       show_col_types = FALSE)
} else {
  message("⚠ Model duyarlılık verisi bulunamadı — ilgili bölümler atlanacak.")
}

# ── GBD ülke katkıları (yeni — ithalat pipeline'dan) ────────────────
COUNTRY_CONTRIB_PATH <- file.path(here("data_processed", SSP_SCENARIO),
                                   "importation_country_contributions.csv")
if (file.exists(COUNTRY_CONTRIB_PATH)) {
  country_contrib <- read_csv(COUNTRY_CONTRIB_PATH, show_col_types = FALSE)
}
```

---

## DEĞİŞİKLİK 3 — Parametre Tablosu [ZORUNLU]

### Satır 175-215: Model Parametreleri tablosu

ρ satırını kaldır veya güncelle. Yeni parametre tablosu:

```r
params_df <- tibble(
  Parametre = c("τ", "IP", "β_vh", "β_hv", "m", "η", "d", "k"),
  Tanım = c(
    "Yerleşme eşiği (majör salgın tanımı)",
    "Enfeksiyöz dönem süresi",
    "Vektör→insan bulaş olasılığı",
    "İnsan→vektör bulaş olasılığı",
    "Sivrisinek/insan oranı (Ross-Macdonald)",
    "Viremik kesir (seyahat sırasında viremik olma olasılığı)",
    "Ortalama maruziyet süresi",
    "İklim elastisite katsayısı"
  ),
  Değer = c("30", "5 gün", "0.30", "0.33", "1.0", "0.25", "3 gün", "0.13 °C⁻¹"),
  Kaynak = c(
    "Operasyonel eşik",
    "Wilder-Smith vd. 2019",
    "Mordecai vd. 2017",
    "Mordecai vd. 2017",
    "Ross-Macdonald standardı",
    "Wilder-Smith vd. 2017, 2019",
    "Lopez vd. 2016; Liebig vd. 2019",
    "Cheng vd. 2023"
  )
)
```

NOT: ρ KALDIRILDI çünkü GBD insidansı zaten under-reporting düzeltmesi
içeriyor. Bunu açıklayan bir dipnot ekle:

> "*GBD 2023 insidans tahminleri, ülke bazlı under-reporting düzeltmesi
> içerdiğinden (Stanaway vd., 2016), ayrı bir ρ (tanı yetersizliği)
> çarpanı kullanılmamıştır.*"

---

## DEĞİŞİKLİK 4 — Yeni Bölüm: GBD Ülke Katkıları [ÖNERİLEN]

Duyarlılık analizi bölümünden (satır 3062) ÖNCE, şu yeni bölümü ekle:

```markdown
## İthalat Baskısının Kaynak Ülke Dağılımı {#sec-ulke-katki}

GBD 2023 ülke bazlı insidans oranları ve TÜİK turist giriş
istatistikleri birleştirilerek, her kaynak ülkenin Türkiye'ye yönelik
dang ithalat riskine katkısı hesaplanmıştır.
```

```r
if (exists("country_contrib")) {
  top_10 <- country_contrib %>%
    slice_max(contribution_pct, n = 10)

  ggplot(top_10, aes(x = reorder(country, contribution_pct),
                     y = contribution_pct)) +
    geom_col(fill = "#E63946") +
    coord_flip() +
    labs(x = NULL, y = "İthalat riskine katkı (%)",
         title = "İlk 10 kaynak ülke — dang ithalat riski",
         subtitle = "GBD 2023 insidans × turist ağırlığı") +
    theme_minimal()
}
```

Bu grafik jüri için çok değerli — "risk nereden geliyor?" sorusuna
somut yanıt veriyor.

---

## DEĞİŞİKLİK 5 — İthalat Duyarlılık Bölümü [ZORUNLU]

### Satır 3062-3100: Eski "η×ρ" referansları → yeni "η + k"

Eski duyarlılık: η×ρ çarpımı (3 düzey) + k (3 düzey)
Yeni duyarlılık: η tek başına (3 düzey) + k (3 düzey)

YENİ BÖLÜM:
```markdown
## İthalat Parametreleri Duyarlılık Analizi {#sec-duyarlilik-ithalat}

İthalat baskısı denklemindeki iki belirsizlik kaynağı — iklim
elastisite katsayısı (k) ve viremik kesir (η) — sistematik olarak
test edilmiştir.
```

```r
if (SENS_IMP_READY) {
  # k duyarlılığı
  k_data <- sens_imp %>% filter(parameter == "k")

  kable(k_data %>%
          select(level, value, year, lambda_total, M_climate) %>%
          pivot_wider(names_from = year, values_from = c(lambda_total, M_climate)),
        caption = "k duyarlılık analizi sonuçları (main senaryo)",
        booktabs = TRUE) %>%
    kable_styling()

  # η duyarlılığı
  eta_data <- sens_imp %>% filter(parameter == "eta")

  kable(eta_data %>%
          select(level, value, year, lambda_total) %>%
          pivot_wider(names_from = year, values_from = lambda_total),
        caption = "η duyarlılık analizi sonuçları (main senaryo)",
        booktabs = TRUE) %>%
    kable_styling()
}
```

---

## DEĞİŞİKLİK 6 — Metrik Sözlüğü Güncellemesi [ZORUNLU]

### Satır 219-260: rho ve pi_global referansları kaldır

ρ yerine:
- pi_weighted: "GBD ülke bazlı insidans ağırlıklı maruz kalma olasılığı"
- rate_pp_year: "GBD ağırlıklı yıllık per-capita insidans oranı"

---

## DEĞİŞİKLİK 7 — Doğrulama Bölümü [KOZMETİK]

### Satır 4300-4421: mc_validation_vs_analytic.R çıktılarını entegre et

Mevcut inline hesaplama yerine, daha önce yazdığımız
`mc_validation_vs_analytic.R` scriptinin çıktılarını yükle:

```r
val_path <- file.path(DIR_DIAG, "mc_validation", "mc_validation_summary.csv")
if (file.exists(val_path)) {
  val_summary <- read_csv(val_path, show_col_types = FALSE)
  kable(val_summary, caption = "MC doğrulama özeti (ilçe bazlı)")
}
```

Ve grafikleri:
```r
val_fig <- file.path(DIR_DIAG, "mc_validation", "fig_mc_vs_analytic.png")
if (file.exists(val_fig)) {
  knitr::include_graphics(val_fig)
}
```

---

## DEĞİŞİKLİK 8 — Özet Bölümü [KOZMETİK]

### Satır 4425-4530: SSP bilgisi ekle

Özet paragrafına SSP senaryosunu belirt:
```r
cat(paste0("Bu sonuçlar ", SSP_LABEL,
           " (middle-of-the-road) senaryosu altında üretilmiştir."))
```

---

## DOKUNULMAYACAK BÖLÜMLER (aynen kalabilir):

- İklim Koşulları (sec-iklim) — iklim verileri SSP-aware olarak
  zaten doğru yoldan yükleniyor
- Yerel Bulaş Kapasitesi (sec-lambda) — λ_local hesabı değişmedi
- Yerleşme Olasılığı (sec-pest) — P_est formülasyonu değişmedi
- Risk Haritası — görsel, veri yapısı aynı
- Uzun Dönemli Risk Trendleri — yearly/horizon veri yapısı aynı
- m Duyarlılık Analizi — sensitivity_ctmc_mc.R yapısı değişmedi
- LHS-PRCC Analizi — bağımsız analiz, etkilenmedi
- CSI Analizi — rainfall_CSI pipeline'ı değişmedi
- Moran's I — bağımsız mekânsal analiz
- Model Şeması — teorik, değişmedi

---

## UYGULAMA ÖNCELİK SIRASI:

1. [ZORUNLU] Dosya yolları (Değişiklik 1)
2. [ZORUNLU] Parametre tablosu - ρ kaldır (Değişiklik 3)
3. [ZORUNLU] Sensitivity veri yükleme (Değişiklik 2)
4. [ZORUNLU] İthalat duyarlılık bölümü (Değişiklik 5)
5. [ÖNERİLEN] GBD ülke katkıları bölümü (Değişiklik 4)
6. [KOZMETİK] Metrik sözlüğü (Değişiklik 6)
7. [KOZMETİK] Doğrulama entegrasyonu (Değişiklik 7)
8. [KOZMETİK] Özet SSP bilgisi (Değişiklik 8)
