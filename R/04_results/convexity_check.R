#!/usr/bin/env Rscript
# =============================================================================
# convexity_check_EIP.R
#
# DOĞRU SORU: Sabit sıcaklıkta, P_est EIP'nin konveks bir fonksiyonu mu?
#
# Gerekçe:
#   Modelinizde Jensen düzeltmesi (make_lambda_local, use_stochastic_EIP=TRUE)
#   SABİT T'de EIP ~ LogNormal(meanlog(T), sdlog(T)) dağılımı üzerinden
#   yapılıyor. Yani Jensen boşluğunun işareti:
#         E[P_est(EIP)]  vs  P_est(E[EIP])
#   EIP dağılımı üzerinden alınır, T dağılımı üzerinden DEĞİL.
#   Dolayısıyla ilgili eğrilik:  d²P_est/dEIP²  (T sabit).
#
# Yöntem:
#   Her ilçe-ay için (temp_c, rh) sabit tutulur; yalnızca EIP oynatılır.
#   λ(EIP) = n·(m·a(T)²·β_vh·β_hv·exp(-μ_v(T,RH)·EIP)) / μ_v(T,RH)
#   R0(EIP) = λ(EIP) at n=1  (spark faz, tek başlangıç enfektesi)
#   P_est(R0) = (1-ρ)/(1-ρ^τ),  ρ = 1/R0   [sonlu-eşikli gambler's ruin]
#   d²P_est/dEIP² merkezi farkla, EIP dağılımının 5.–95. persentil
#   aralığında (yani Jensen'in fiilen örneklediği bölgede) taranır.
#
# Çıktı:
#   - Her ilçe-ay-yıl için baskın eğrilik (EIP dağılım kütlesi ağırlıklı)
#   - Jensen boşluğunun İŞARETİ ile eğriliğin tutarlılık kontrolü
#   - İşaret haritası + gerçek dağılım aralıklarının bindirilmesi
#
# Okan Derin — dengue CTMC tezi
# Çalıştırma: R/03_models altından source("parameter_functions.R") erişimiyle
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(readr)
})

# ---- 0. YOLLAR (kendi projenize göre; .Rproj kökünden) ---------------------
source("R/03_models/parameter_functions.R")

TRAIT_PATH <- "data_processed/trait_params.csv"   # tek tür ise
# Tür bazlı ayrı tablolar kullanıyorsanız aşağıda species'e göre seçilecek:
TRAIT_PATH_AEGYPTI    <- "data_processed/trait_params_aegypti.csv"
TRAIT_PATH_ALBOPICTUS <- "data_processed/trait_params_albopictus.csv"

TAU <- 30

# ---- 1. SABİT TRANSMİSYON PARAMETRELERİ ------------------------------------
# ctmc_spark.R içinde kullandığınız değerlerle BİREBİR eşleştirin.
# Aşağıdakiler yer tutucu; kendi değerlerinizi girin ya da tablodan okuyun.
BETA_VH <- 0.5     # vector->human
BETA_HV <- 0.5     # human->vector
M_RATIO <- 1.0     # m (mosquito/human) — CSI/mevsimsel modülasyon öncesi taban

# ---- 2. P_est (sonlu-tau gambler's ruin) -----------------------------------
Pest_finite <- function(R0, tau = TAU) {
  rho <- 1 / R0
  out <- ifelse(abs(R0 - 1) < 1e-9, 1/tau, (1 - rho) / (1 - rho^tau))
  out[R0 <= 0 | !is.finite(R0)] <- 0
  pmin(pmax(out, 0), 1)
}

# ---- 3. R0(EIP) sabit T,RH'de (gerçek trait fonksiyonlarıyla) --------------
# EIP'yi burada BAĞIMSIZ değişken olarak veriyoruz; a(T), mu_v(T,RH) T'den gelir.
R0_of_EIP_at_T <- function(EIP, T, RH, tab, beta_vh = BETA_VH,
                           beta_hv = BETA_HV, m = M_RATIO) {
  aT  <- a_of_T(T, tab)
  muv <- mu_v_of_TRH(T, RH, tab)
  surv <- exp(-muv * EIP)
  (m * aT^2 * beta_vh * beta_hv * surv) / muv   # n=1
}

Pest_of_EIP_at_T <- function(EIP, T, RH, tab, ...) {
  Pest_finite(R0_of_EIP_at_T(EIP, T, RH, tab, ...))
}

# İkinci türev (merkezi fark), sabit T,RH
d2_Pest_dEIP2 <- function(EIP, T, RH, tab, h = 1e-3, ...) {
  fp <- Pest_of_EIP_at_T(EIP + h, T, RH, tab, ...)
  f0 <- Pest_of_EIP_at_T(EIP,     T, RH, tab, ...)
  fm <- Pest_of_EIP_at_T(EIP - h, T, RH, tab, ...)
  (fp - 2*f0 + fm) / h^2
}

# ---- 4. GERÇEK VERİYİ YÜKLE (üç SSP birleşik) ------------------------------
load_scenario <- function(ssp) {
  readRDS(sprintf("outputs/%s/model_results/ctmc_spark_monthly_2025_2075.rds", ssp)) %>%
    mutate(ssp = ssp)
}
dat <- bind_rows(lapply(c("ssp126","ssp245","ssp585"), load_scenario))

# Tür bazlı trait tablosu seçici (species kolonu var)
get_tab <- local({
  cache <- list()
  function(species) {
    key <- tolower(species)
    if (is.null(cache[[key]])) {
      p <- if (key == "aegypti") TRAIT_PATH_AEGYPTI else TRAIT_PATH_ALBOPICTUS
      if (!file.exists(p)) p <- TRAIT_PATH   # tek tablo fallback
      cache[[key]] <<- read_trait_params(p)
    }
    cache[[key]]
  }
})

# ---- 5. HER SATIR İÇİN: EIP dağılım aralığı + eğrilik ----------------------
# Jensen'in fiilen örneklediği aralık: LogNormal 5.–95. persentil.
# Eğriliği o aralıkta değerlendiriyoruz (uçlarda değil).
analyze_row <- function(temp_c, rh, species) {
  tab <- get_tab(species)
  
  # EIP dağılım parametreleri (sizin kalibre fonksiyonunuz)
  pars <- get_eip_logn_params(temp_c, tab, force_recalibrate = FALSE)
  if (!is.finite(pars$meanlog)) {
    return(tibble(eip_mean = NA_real_, eip_q05 = NA_real_, eip_q95 = NA_real_,
                  d2_at_mean = NA_real_, frac_convex = NA_real_,
                  egrilik = "termal_pencere_disi"))
  }
  
  eip_mean <- exp(pars$meanlog + 0.5 * pars$sdlog^2)
  eip_q05  <- qlnorm(0.05, pars$meanlog, pars$sdlog)
  eip_q95  <- qlnorm(0.95, pars$meanlog, pars$sdlog)
  
  # Dağılım kütlesinin yoğun olduğu aralıkta d² işaretini tara
  eip_grid <- seq(eip_q05, eip_q95, length.out = 41)
  d2_grid  <- vapply(eip_grid, function(e)
    d2_Pest_dEIP2(e, temp_c, rh, tab), numeric(1))
  
  # Dağılım yoğunluğuyla ağırlıklandır (Jensen boşluğuna asıl katkı burdan)
  w <- dlnorm(eip_grid, pars$meanlog, pars$sdlog)
  
  # NA/non-finite d² noktalarını (termal pencere sınırı, R0=0 vb.) DIŞLA
  ok <- is.finite(d2_grid) & is.finite(w)
  if (!any(ok) || sum(w[ok]) <= 0) {
    return(tibble(eip_mean = eip_mean, eip_q05 = eip_q05, eip_q95 = eip_q95,
                  d2_at_mean = NA_real_, frac_convex = NA_real_,
                  egrilik = "belirsiz"))
  }
  w  <- w[ok] / sum(w[ok])
  d2 <- d2_grid[ok]
  frac_convex <- sum(w[d2 > 0])
  
  d2_at_mean <- d2_Pest_dEIP2(eip_mean, temp_c, rh, tab)
  
  egrilik <- if (!is.finite(frac_convex)) "belirsiz"
  else if (frac_convex > 0.99) "Konveks"
  else if (frac_convex < 0.01) "Konkav"
  else "Karisik"
  
  tibble(eip_mean, eip_q05, eip_q95, d2_at_mean, frac_convex, egrilik)
}

res <- dat %>%
  mutate(row = pmap(list(temp_c, rh, species), analyze_row)) %>%
  unnest(row)

# ---- 6. ÖZET ----------------------------------------------------------------
cat("\n=== Eğrilik dağılımı (EIP dağılım-ağırlıklı, sabit T) ===\n")
res %>% count(egrilik) %>% mutate(pct = round(100*n/sum(n),1)) %>% print()

# "belirsiz" ve "termal_pencere_disi" tipik olarak kış/soğuk aylar:
# EIP sonsuz veya a(T)=0, yani zaten p_establishment≈0. Bunları
# aktif (p_establishment>0) satırlar arasında olup olmadığını kontrol et:
cat("\n=== Analiz-dışı satırların aktiflik kontrolü ===\n")
res %>%
  filter(egrilik %in% c("belirsiz","termal_pencere_disi")) %>%
  summarise(n_toplam = n(),
            n_aktif  = sum(p_establishment > 0, na.rm = TRUE)) %>%
  print()
cat("(n_aktif=0 beklenir: analiz dışı kalanlar zaten sıfır-risk aylardır.)\n")

cat("\n=== SSP × ilçe bazında konveks kütle oranı (ortalama) ===\n")
res %>%
  filter(!is.na(frac_convex)) %>%
  group_by(ssp, district_name.x) %>%
  summarise(mean_frac_convex = round(mean(frac_convex),3),
            n_active = sum(p_establishment > 0), .groups = "drop") %>%
  arrange(ssp, desc(mean_frac_convex)) %>%
  print(n = Inf)

# ---- 7. JENSEN İŞARET TUTARLILIK KONTROLÜ ----------------------------------
# validation_stage2_jensen.csv ile çapraz doğrulama:
#   frac_convex ~ 1  => E[P_est] >= P_est(E)  => jensen_gap >= 0 beklenir
jf <- tryCatch(
  bind_rows(lapply(c("ssp126","ssp245","ssp585"), function(s)
    read_csv(sprintf("outputs/%s/validation/validation_stage2_jensen.csv", s),
             show_col_types = FALSE) %>% mutate(ssp = s))),
  error = function(e) { message("Jensen CSV okunamadı: ", e$message); NULL })

if (!is.null(jf)) {
  cat("\n=== Jensen boşluğu işareti vs eğrilik (beklenen: konveks=pozitif gap) ===\n")
  cat("validation_stage2_jensen.csv kolonları:\n  ",
      paste(names(jf), collapse = ", "), "\n")
  cat("Bu tabloyu res ile district/year/month üzerinden join edip\n",
      "sign(jensen_gap) == (frac_convex>0.5) tutarlılığını kontrol edin.\n")
}

# ---- 8. İŞARET HARİTASI (T ekseni) + gerçek EIP dağılımları -----------------
# Sabit-T eğrilik, T'ye göre nasıl değişiyor? Her tür için ayrı.
plot_species <- "aegypti"
tab_p <- get_tab(plot_species)
Tgrid <- seq(15, 38, length.out = 200)
map_df <- map_dfr(Tgrid, function(T) {
  pars <- get_eip_logn_params(T, tab_p)
  if (!is.finite(pars$meanlog)) return(NULL)
  eip_mean <- exp(pars$meanlog + 0.5*pars$sdlog^2)
  d2 <- d2_Pest_dEIP2(eip_mean, T, 70, tab_p)  # tipik RH=70
  if (!is.finite(d2)) return(NULL)
  tibble(T, eip_mean, d2, egrilik = ifelse(d2 >= 0, "Konveks", "Konkav"))
})

p_map <- ggplot(map_df, aes(T, eip_mean, color = egrilik)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Konveks"="#2C7FB8","Konkav"="#D7301F")) +
  labs(title = sprintf("d²P_est/dEIP² işareti — %s (RH=70)", plot_species),
       subtitle = "EIP ortalamasında değerlendirildi; sabit T",
       x = "Sıcaklık (°C)", y = "Ortalama EIP (gün)", color = "Eğrilik") +
  theme_minimal(base_size = 12)
ggsave("outputs/convexity_EIP_map.png", p_map, width = 8, height = 5, dpi = 150)

# ---- 9. TEMSİLİ EĞRİ: aktif bir ilçe-ay için P_est(EIP) --------------------
active <- dat %>% filter(p_establishment > 0.01) %>% 
  arrange(desc(p_establishment)) %>% slice(1)
if (nrow(active) == 1) {
  tab_a <- get_tab(active$species)
  pars  <- get_eip_logn_params(active$temp_c, tab_a)
  eE    <- exp(pars$meanlog + 0.5*pars$sdlog^2)
  eip_seq <- seq(max(0.5, qlnorm(0.01,pars$meanlog,pars$sdlog)),
                 qlnorm(0.99,pars$meanlog,pars$sdlog), length.out = 300)
  curve_df <- tibble(
    EIP  = eip_seq,
    Pest = Pest_of_EIP_at_T(eip_seq, active$temp_c, active$rh, tab_a),
    dens = dlnorm(eip_seq, pars$meanlog, pars$sdlog)
  )
  p_curve <- ggplot(curve_df, aes(EIP, Pest)) +
    geom_line(color = "#2C7FB8", linewidth = 0.9) +
    geom_vline(xintercept = eE, linetype = "dashed") +
    geom_area(aes(y = dens/max(dens)*max(Pest)*0.4), fill = "grey70", alpha = 0.4) +
    labs(title = sprintf("Temsili aktif ilçe-ay: %s, %.1f°C",
                         active$district_name.x, active$temp_c),
         subtitle = "Kesikli çizgi = E[EIP]; gri = EIP dağılımı (ölçekli)",
         x = "EIP (gün)", y = expression(P[est])) +
    theme_minimal(base_size = 12)
  ggsave("outputs/convexity_EIP_curve.png", p_curve, width = 8, height = 5, dpi = 150)
}

cat("\nGrafikler: outputs/convexity_EIP_map.png, outputs/convexity_EIP_curve.png\n")
cat("Eğrilik özet tablosu 'res' nesnesinde. Rapor için:\n")
cat("  readr::write_csv(res, 'outputs/convexity_EIP_summary.csv')\n")