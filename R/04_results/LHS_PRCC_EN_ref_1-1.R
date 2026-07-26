## =========================================================
## LHS–PRCC duyarlılık grafiği  —  x-ekseni [-1, +1]
## 0 referans çizgisi + işaret-bazlı renk (poz/neg)
## Panel düzeni orijinaldeki gibi: 2x2, sağ-alt boş
## Gereksinim: ggplot2, patchwork  (forcats opsiyonel)
## =========================================================

library(ggplot2)
library(patchwork)

## ---------------------------------------------------------
## Beklenen veri: her satır bir (senaryo, parametre) PRCC sonucu
##   scenario  : "SSP1-2.6" / "SSP2-4.5" / "SSP5-8.5"
##   parameter : "T_C","m","beta_vh","beta_hv","RH","ip_days"
##   prcc      : nokta PRCC tahmini (∈ [-1, 1])
##   conf_lo,conf_hi : bootstrap %95 GA
## Alttaki ÖRNEK, grafiğinizdeki değerlere yakındır —
## kendi prcc_results tablonuzla değiştirin.
## ---------------------------------------------------------
prcc_results <- data.frame(
  scenario  = rep(c("SSP1-2.6","SSP2-4.5","SSP5-8.5"), each = 6),
  parameter = rep(c("T_C","m","beta_vh","beta_hv","RH","ip_days"), 3),
  prcc      = c(0.90,0.67,0.54,0.51,0.51,0.42,
                0.88,0.68,0.55,0.55,0.50,0.42,
                0.90,0.65,0.53,0.52,0.51,0.40),
  conf_lo   = c(0.87,0.63,0.49,0.46,0.46,0.37,
                0.85,0.64,0.50,0.50,0.45,0.37,
                0.87,0.61,0.48,0.47,0.46,0.35),
  conf_hi   = c(0.92,0.71,0.59,0.56,0.56,0.47,
                0.90,0.72,0.60,0.60,0.55,0.47,
                0.92,0.69,0.58,0.57,0.56,0.45),
  stringsAsFactors = FALSE
)

## Parametreleri sabit sırada tut (T_C en üstte) + okunur etiketler ----
disp <- c(ip_days = "Bulaşıcı dönem (ip_days)",
          RH      = "Bağıl nem (RH)",
          beta_hv = "\u03b2_hv",
          beta_vh = "\u03b2_vh",
          m       = "Sivrisinek:insan (m)",
          T_C     = "Sıcaklık (T_C)")
lev <- c("ip_days","RH","beta_hv","beta_vh","m","T_C")   # alt -> üst
prcc_results$parameter <- factor(disp[prcc_results$parameter],
                                 levels = unname(disp[lev]))

## ---------------------------------------------------------
## Tek senaryo için PRCC grafiği (x-ekseni [-1, +1])
## ---------------------------------------------------------
plot_prcc <- function(df, title = NULL, subtitle = NULL,
                      pos_col = "#E24A4A", neg_col = "#3B6FB6") {
  ggplot(df, aes(x = prcc, y = parameter, fill = prcc >= 0)) +
    geom_col(width = 0.70) +
    geom_errorbarh(aes(xmin = conf_lo, xmax = conf_hi),
                   height = 0.25, linewidth = 0.4, colour = "grey20") +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +   # ← referans
    scale_fill_manual(values = c(`TRUE` = pos_col, `FALSE` = neg_col),
                      guide = "none") +
    scale_x_continuous(limits = c(-1, 1),                              # ← TAM ARALIK
                       breaks = seq(-1, 1, 0.5),
                       expand = expansion(mult = 0.02)) +
    labs(x = "PRCC", y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, colour = "grey45"),
      axis.text.y   = element_text(size = 9)
    )
}

## Panel alt başlıkları (T/RH aralıkları senaryoya özel) ----
subs <- c(
  "SSP1-2.6" = "n = 2000 | T: 10,4\u201329,6\u00b0C | RH: 21,0\u201396,2% (active-month bootstrap)",
  "SSP2-4.5" = "n = 2000 | T: 10,4\u201330,0\u00b0C | RH: 10,5\u201396,4% (active-month bootstrap)",
  "SSP5-8.5" = "n = 2000 | T: 10,4\u201332,1\u00b0C | RH: 18,1\u201397,1% (active-month bootstrap)"
)

## Üç paneli üret ----
plots <- lapply(names(subs), function(s) {
  plot_prcc(prcc_results[prcc_results$scenario == s, ],
            title    = paste(s, "LHS\u2013PRCC Sensitivity Analysis"),
            subtitle = subs[[s]])
})

## 2x2 düzen, sağ-alt boş (orijinaldeki gibi) ----
fig_prcc <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plot_spacer())

## Kaydet ----
ggsave("fig_prcc_pm1.png", fig_prcc, width = 11, height = 7, dpi = 300)

## =========================================================
## ALTERNATİF: patchwork yerine facet kullanıyorsanız,
## tek gereken değişiklik x-ekseni + 0 çizgisidir:
##
##   ... + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +
##         scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
##         facet_wrap(~ scenario)
##
## (facet ile panel-başına farklı alt başlık veremezsiniz;
##  o yüzden T/RH aralıkları için patchwork tercih edilmiştir.)
## =========================================================