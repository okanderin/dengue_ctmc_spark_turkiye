rm(list = ls())

# Adım 1: Her SSP için grafik ve tabloları üret (~2-3 dk/SSP)
for (s in c("ssp126", "ssp245", "ssp585")) {
  Sys.setenv(SSP_SCENARIO =s)
  source("R/04_results/01_generate_ssp_outputs.R", local = TRUE)
}


# Adım 2: Karşılaştırmalı raporu knit et (~30 sn)
rmarkdown::render("R/04_results/bulgular.Rmd")


# Yeni Rmd'nin ürettiği cross-SSP çıktılar:

## Yıllık risk trendleri — 3 SSP yan yana, ilçe bazlı facet, güven aralıkları ile. 
## 50 yıllık ufuk tablosu — her ilçe × 3 SSP matrisinde P_ufuk, GA, Λ_ithalat. 
## Isınma etkisi — erken vs geç dönem ΔT ve kat değişim karşılaştırması. 
## Kaynak ülke analizi — GBD ülke katkıları grafiği (pipeline'ın en güçlü yeni çıktısı). 
## İthalat duyarlılığı — k ve η tablosu.

# Eski Rmd'den korunan ve 01_generate_ssp_outputs.R'ye taşınan analizler:

## İklim profilleri (sıcaklık/nem mevsimselliği), 
## λ_local ısı haritası, 
## P_est mevsimselliği, 
## yıllık risk trendi, 
## ufuk tablosu,
## model doğrulaması, 
## tornado duyarlılık grafiği, 
## ülke katkıları. 

## LHS-PRCC analizi, 
## CSI ısı haritası ve trend, 
## Moran's I mekânsal otokorelasyon, 
## dekadal risk progresyonu, 
## bulaş sezonu uzunluğu değişimi, 
## çoklu regresyon (T + RH + Λ_import → risk).
