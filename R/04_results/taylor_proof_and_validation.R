## =========================================================
## Taylor Açılımı ile P_ufuk Yakınsama Kanıtı
## Tez Metodoloji Bölümüne Eklenecek Metin + Doğrulama Scripti
## =========================================================

## =============================================================
## BÖLÜM 1: TEZ METNİ (Türkçe, LaTeX formatlı)
## =============================================================

# ---- Önerilen bölüm başlığı: "3.X.X Kümülatif Risk ve Hazard İntegrali Yakınsaması"

# --- TEZ METNİ BAŞLANGIÇ ---
#
# Modelin birincil çıktısı olan P_ufuk, projeksiyon dönemi boyunca en az bir
# ayda otokton yerleşmenin gerçekleşme olasılığını ifade etmektedir. Bu değer,
# aylık salgın olasılıklarının tümleyen çarpımı olarak hesaplanmaktadır:
#
#   P_ufuk = 1 - ∏_{i=1}^{N} (1 - p_i)                           (12)
#
# burada p_i = q_ithal,i × P_est,i aylık salgın olasılığı, N = 612
# (51 yıl × 12 ay) toplam ay sayısıdır. Bu formülasyon, aylar arası
# bağımsızlık varsayımı altında geçerlidir.
#
# Denklem (12), sürekli zamanlı kümülatif hazard integrali ile
# ilişkilendirilebilir. Logaritmik dönüşüm uygulandığında:
#
#   ln(1 - P_ufuk) = ∑_{i=1}^{N} ln(1 - p_i)                    (12a)
#
# Taylor serisinin birinci mertebe açılımına göre, |p_i| < 1 koşulunda:
#
#   ln(1 - p_i) = -p_i - p_i²/2 - p_i³/3 - ...                  (12b)
#
# p_i ≪ 1 olduğunda yüksek mertebe terimleri ihmal edilebilir ve:
#
#   ln(1 - p_i) ≈ -p_i                                            (12c)
#
# Bu yaklaşımın hata sınırı ikinci mertebe terimden gelir:
#
#   |R_i| = |ln(1 - p_i) + p_i| = p_i²/2 + O(p_i³) ≤ p_i²/2    (12d)
#
# Yaklaşım (12c) Denklem (12a)'ya uygulandığında:
#
#   ln(1 - P_ufuk) ≈ -∑_{i=1}^{N} p_i = -Λ                      (12e)
#
# Dolayısıyla:
#
#   P_ufuk ≈ 1 - exp(-Λ)                                          (12f)
#
# Bu, kümülatif hazard Λ = ∑ p_i ile tanımlanan Poisson süreç
# yaklaşımının sonucudur ve sürekli zamanlı hazard integralinin
# (P = 1 - exp(-∫h(t)dt)) ayrık zaman karşılığıdır.
#
# Yaklaşımın geçerlilik koşulu: toplam hata sınırı
#
#   |R_toplam| ≤ ∑_{i=1}^{N} p_i²/2 = (1/2) ∑ p_i²              (12g)
#
# olmak üzere, bu değerin ihmal edilebilir düzeyde kalması için
# p_i değerlerinin küçük olması yeterlidir. Çalışmamızda beş
# sentinel ilçe ve üç SSP senaryosu için bu koşulun sağlanma
# durumu Tablo XX'de gösterilmektedir.
#
# --- TEZ METNİ BİTİŞ ---


## =============================================================
## BÖLÜM 2: DOĞRULAMA SCRİPTİ
## Kesin formül vs Poisson yaklaşımı karşılaştırması
## =============================================================

# Bu script, bulgular.Rmd'ye chunk olarak eklenecektir.
# Ön koşul: all_monthly listesi yüklenmiş olmalı.

# ```{r taylor-validation, fig.cap="Kesin tümleyen çarpımı ile Poisson hazard yaklaşımının karşılaştırması."}

library(tidyverse)
library(scales)

Sys.setenv(SSP_SCENARIO = "ssp585")

## ---- 1) Her ilçe-SSP için kesin vs Poisson karşılaştırması ----
taylor_results <- list()

for (s in SSP_LIST) {
  dat_s <- all_monthly[[s]]
  if (is.null(dat_s) || nrow(dat_s) == 0) next
  
  taylor_results[[s]] <- dat_s %>%
    mutate(
      # p_ay = q_import × P_est (aylık salgın olasılığı)
      p_ay = pmin(pmax(p_month_major_mean, 0), 1)
    ) %>%
    group_by(district_id, district_label) %>%
    summarise(
      # Kesin formül: P = 1 - ∏(1 - p_i)
      P_exact = 1 - prod(1 - p_ay, na.rm = TRUE),
      
      # Poisson yaklaşımı: P ≈ 1 - exp(-Σp_i)
      Lambda_cum = sum(p_ay, na.rm = TRUE),
      P_poisson = 1 - exp(-Lambda_cum),
      
      # Taylor hata diagnostikleri
      max_p_ay = max(p_ay, na.rm = TRUE),
      mean_p_ay = mean(p_ay, na.rm = TRUE),
      n_months = n(),
      
      # Aktif ay sayısı (p_ay > 0)
      n_active = sum(p_ay > 0),
      
      # Toplam Taylor hata sınırı: (1/2) Σ p_i²
      R_total = 0.5 * sum(p_ay^2, na.rm = TRUE),
      
      # p > 0.01 olan ay sayısı (Taylor yaklaşımı zayıflayan aylar)
      n_large_p = sum(p_ay > 0.01, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(
      # Kesin vs Poisson farkı
      abs_diff = abs(P_exact - P_poisson),
      rel_diff_pct = ifelse(P_exact > 0,
                            abs_diff / P_exact * 100,
                            0),
      ssp = s,
      ssp_label = SSP_LABELS[s]
    )
}

taylor_all <- bind_rows(taylor_results)

## ---- 2) Özet tablo: Taylor yakınsama diagnostiği ----
taylor_summary <- taylor_all %>%
  select(
    Senaryo    = ssp_label,
    İlçe       = district_label,
    `P_ufuk (kesin)`    = P_exact,
    `P_ufuk (Poisson)`  = P_poisson,
    `Λ (küm. hazard)`   = Lambda_cum,
    `maks p_ay`          = max_p_ay,
    `Σp²/2 (hata sınırı)` = R_total,
    `Bağıl fark (%)`     = rel_diff_pct,
    `p > 0.01 ay sayısı` = n_large_p
  ) %>%
  arrange(Senaryo, desc(`P_ufuk (kesin)`))

# Tablo formatı
taylor_display <- taylor_summary %>%
  mutate(
    across(c(`P_ufuk (kesin)`, `P_ufuk (Poisson)`),
           ~formatC(.x, format = "e", digits = 3)),
    `Λ (küm. hazard)` = round(`Λ (küm. hazard)`, 4),
    `maks p_ay` = formatC(`maks p_ay`, format = "e", digits = 2),
    `Σp²/2 (hata sınırı)` = formatC(`Σp²/2 (hata sınırı)`, format = "e", digits = 2),
    `Bağıl fark (%)` = round(`Bağıl fark (%)`, 4)
  )

kable(taylor_display,
      caption = paste0(
        "Taylor yakınsama doğrulaması: Kesin tümleyen çarpımı (Denklem 12) ile ",
        "Poisson kümülatif hazard yaklaşımı (Denklem 12f) karşılaştırması. ",
        "Σp²/2, Taylor açılımının ikinci mertebe hata sınırını göstermektedir."
      ),
      booktabs = TRUE,
      escape = FALSE) %>%
  kable_styling(font_size = 10, full_width = FALSE)


## ---- 3) Kartal detay analizi: hangi aylar Taylor'ı ihlal ediyor? ----
kartal_detail <- NULL

for (s in SSP_LIST) {
  dat_s <- all_monthly[[s]]
  if (is.null(dat_s) || nrow(dat_s) == 0) next
  
  kartal_s <- dat_s %>%
    filter(grepl("Kartal|TUR.40.25", district_id) |
             grepl("Kartal|TUR.40.25", district_label)) %>%
    mutate(
      p_ay = pmin(pmax(p_month_major_mean, 0), 1),
      taylor_error = p_ay^2 / 2,
      taylor_valid = p_ay < 0.01,  # %1 eşiği
      ssp = s,
      ssp_label = SSP_LABELS[s]
    ) %>%
    filter(p_ay > 0) %>%
    select(ssp_label, year, month, p_ay, taylor_error, taylor_valid,
           q_import_month, p_establishment_mean)
  
  kartal_detail <- bind_rows(kartal_detail, kartal_s)
}

# Kartal'da Taylor yaklaşımının zayıfladığı aylar
kartal_violations <- kartal_detail %>%
  filter(!taylor_valid) %>%
  arrange(ssp_label, desc(p_ay))

cat("\n=== Kartal: Taylor yaklaşımının zayıfladığı aylar (p_ay > 0.01) ===\n")
cat("Toplam ihlal sayısı:", nrow(kartal_violations), "/",
    nrow(kartal_detail), "aktif ay\n\n")

if (nrow(kartal_violations) > 0) {
  kartal_violation_summary <- kartal_violations %>%
    group_by(ssp_label) %>%
    summarise(
      n_ihlal = n(),
      max_p_ay = max(p_ay),
      mean_p_ay_ihlal = mean(p_ay),
      toplam_hata = sum(taylor_error),
      tipik_aylar = paste(unique(month), collapse = ", "),
      .groups = "drop"
    )
  print(kartal_violation_summary)
}


## ---- 4) Görselleştirme: Kesin vs Poisson scatter ----
fig_taylor <- ggplot(taylor_all,
                     aes(x = P_exact, y = P_poisson,
                         colour = district_label,
                         shape = ssp_label)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_x_log10(labels = label_scientific()) +
  scale_y_log10(labels = label_scientific()) +
  scale_colour_manual(values = COL_DISTRICT, name = "İlçe") +
  scale_shape_manual(values = c(16, 17, 15), name = "Senaryo") +
  labs(x = expression(P[ufuk]^{kesin}~"(tümleyen çarpımı)"),
       y = expression(P[ufuk]^{Poisson}~"(hazard yaklaşımı)")) +
  theme_thesis() +
  # Kartal etiketlerini ekle (potansiyel sapma noktaları)
  annotate("text", x = max(taylor_all$P_exact) * 0.3,
           y = max(taylor_all$P_poisson) * 0.1,
           label = "Birim eğim: kesin = Poisson",
           colour = "red", size = 3, fontface = "italic")

print(fig_taylor)

ggsave(
  filename = file.path(DIR_OUTPUT_CROSS, "fig_taylor_validation.png"),
  plot = fig_taylor,
  width = 8, height = 6, dpi = 300
)


## ---- 5) Kartal aylik p_ay profili (hangi aylarda büyük?) ----
kartal_monthly_profile <- kartal_detail %>%
  group_by(ssp_label, month) %>%
  summarise(
    p_ay_mean = mean(p_ay, na.rm = TRUE),
    p_ay_max  = max(p_ay, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(month_name = factor(AY_TR[month], levels = AY_TR))

fig_kartal_profile <- ggplot(kartal_monthly_profile,
                              aes(x = month_name, y = p_ay_max,
                                  fill = ssp_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0.01, linetype = "dashed",
             colour = "red", linewidth = 0.8) +
  annotate("text", x = 1, y = 0.013, label = "Taylor eşiği (p = 0.01)",
           colour = "red", size = 3, hjust = 0) +
  scale_fill_manual(values = COL_SSP, name = "Senaryo") +
  scale_y_continuous(labels = label_scientific()) +
  labs(x = NULL, y = expression(p[ay]~"(maks)"),
       title = "Kartal — aylık salgın olasılığı profili") +
  theme_thesis()

print(fig_kartal_profile)

ggsave(
  filename = file.path(DIR_OUTPUT_CROSS, "fig_kartal_taylor_profile.png"),
  plot = fig_kartal_profile,
  width = 8, height = 5, dpi = 300
)

# ```


## =============================================================
## BÖLÜM 3: KARTAL'IN TAYLOR YAKLAŞIMINI NEDEN İHLAL ETTİĞİ
##          (Tez tartışma bölümüne eklenecek açıklama)
## =============================================================

# --- TEZ METNİ ---
#
# Kartal ilçesi, Taylor yaklaşımının kısmen zayıfladığı tek sentinel
# ilçedir. Bu durumun nedeni, Kartal'ın iki koşulu aynı anda
# sağlamasıdır:
#
# (i)  Yüksek ithalat baskısı: Kartal, İstanbul'un uluslararası
#      havalimanlarına yakınlığı nedeniyle yüksek q_ithal değerlerine
#      sahiptir (aylık ithal vaka olasılığı).
#
# (ii) Yüksek yerel bulaş kapasitesi: Kartal'ın yazın (Temmuz–Eylül)
#      sıcaklık ve nem koşulları Ae. albopictus için termal optimuma
#      yakındır, bu da P_est değerini ~10⁻¹ mertebesine çıkarmaktadır.
#
# Bu iki faktörün çarpımı olan p_ay = q_ithal × P_est, yaz aylarında
# 10⁻² eşiğini aşabilmektedir. Bu düzeyde Taylor açılımının birinci
# mertebe yaklaşımı (ln(1-p) ≈ -p) artık yeterli doğrulukta değildir;
# ikinci mertebe terimi (p²/2) anlamlı hale gelmektedir.
#
# Ancak bu durum modelin geçerliliğini zayıflatmamaktadır, çünkü:
# (a) P_ufuk hesabında kesin tümleyen çarpımı (Denklem 12) kullanılmakta
#     olup Poisson yaklaşımı bir hesaplama aracı olarak değil, yalnızca
#     kavramsal çerçeve ve doğrulama amaçlı sunulmaktadır.
# (b) Kartal için bile kesin ve Poisson tahminleri arasındaki bağıl fark
#     %X düzeyinde kalmaktadır (Tablo XX).
# (c) Diğer dört ilçede tüm aylarda p_ay < 10⁻² olup Taylor yaklaşımı
#     tam olarak geçerlidir.
#
# --- TEZ METNİ BİTİŞ ---
