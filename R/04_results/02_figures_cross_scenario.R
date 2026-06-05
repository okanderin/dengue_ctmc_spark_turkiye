
library(magick)
library(patchwork)
library(ggplot2)
library(grid)

# SSP listesini ve çıktı klasörünü tanımlayın
SSP_LIST <- c("ssp126", "ssp245", "ssp585") # Klasör isimlerinizle aynı olmalı

# 1. Resimleri oku ve ggplot nesnesine çevir
plot_list_temp_rh <- lapply(SSP_LIST, function(s) {
  path <- here("outputs", s, "figures", "fig_temp_rh.png")
  if (file.exists(path)) {
    img <- image_read(path)
    return(rasterGrob(img, interpolate = TRUE))
  }
  return(NULL)
})

# NULL olanları temizle
plot_list_temp_rh <- plot_list_temp_rh[!sapply(plot_list_temp_rh, is.null)]

# 2. Patchwork ile birleştir (Örn: Yan yana 3 grafik)
combined_plot_temp_rh_profile <- wrap_plots(plot_list_temp_rh) + 
  plot_layout(ncol = 2)

# 3. 300 DPI olarak kaydet
ggsave(
  filename = here(DIR_OUTPUT_CROSS, "combined_plot_temp_rh_profile.png"),
  plot = combined_plot_temp_rh_profile,
  dpi = 300,
  width = 15, # İnç cinsinden genişlik
  height = 10,  # İnç cinsinden yükseklik
  units = "in"
)




# 2. CSI heatmap Resimleri oku ve ggplot nesnesine çevir
plot_list_csi_heat <- lapply(SSP_LIST, function(s) {
  path <- here("outputs", s, "figures", "fig_csi_heat.png")
  if (file.exists(path)) {
    img <- image_read(path)
    return(rasterGrob(img, interpolate = TRUE))
  }
  return(NULL)
})

# NULL olanları temizle
plot_list_csi_heat <- plot_list_csi_heat[!sapply(plot_list_csi_heat, is.null)]

# 2. Patchwork ile birleştir (Örn: Yan yana 3 grafik)
combined_plot_csi_heat <- wrap_plots(plot_list_csi_heat) + 
  plot_layout(ncol = 2)

# 3. 300 DPI olarak kaydet
ggsave(
  filename = here(DIR_OUTPUT_CROSS, "combined_plot_csi_heat.png"),
  plot = combined_plot_csi_heat,
  dpi = 300,
  width = 15, # İnç cinsinden genişlik
  height = 10,  # İnç cinsinden yükseklik
  units = "in"
)


library(magick)
library(grid)
library(here)


# 3. CSI yıllık tend Resimleri oku ve ggplot nesnesine çevir
plot_list_csi_trend <- lapply(SSP_LIST, function(s) {
  path <- here("outputs", s, "figures", "fig_csi_trend.png")
  if (file.exists(path)) {
    img <- image_read(path)
    return(rasterGrob(img, interpolate = TRUE))
  }
  return(NULL)
})

# NULL olanları temizle
plot_list_csi_trend <- plot_list_csi_trend[!sapply(plot_list_csi_trend, is.null)]

# 2. Patchwork ile birleştir (Örn: Yan yana 3 grafik)
combined_plot_csi_trend <- wrap_plots(plot_list_csi_trend) + 
  plot_layout(ncol = 2)

# 3. 300 DPI olarak kaydet
ggsave(
  filename = here(DIR_OUTPUT_CROSS, "combined_plot_csi_trend.png"),
  plot = combined_plot_csi_trend,
  dpi = 300,
  width = 15, # İnç cinsinden genişlik
  height = 10,  # İnç cinsinden yükseklik
  units = "in"
)



# 4.Aylık risk Resimleri oku ve ggplot nesnesine çevir
plot_list_csi_trend <- lapply(SSP_LIST, function(s) {
  path <- here("outputs", s, "figures", "fig_csi_trend.png")
  if (file.exists(path)) {
    img <- image_read(path)
    return(rasterGrob(img, interpolate = TRUE))
  }
  return(NULL)
})

# NULL olanları temizle
plot_list_csi_trend <- plot_list_csi_trend[!sapply(plot_list_csi_trend, is.null)]

# 2. Patchwork ile birleştir (Örn: Yan yana 3 grafik)
combined_plot_csi_trend <- wrap_plots(plot_list_csi_trend) + 
  plot_layout(ncol = 2)

# 3. 300 DPI olarak kaydet
ggsave(
  filename = here(DIR_OUTPUT_CROSS, "combined_plot_csi_trend.png"),
  plot = combined_plot_csi_trend,
  dpi = 300,
  width = 15, # İnç cinsinden genişlik
  height = 10,  # İnç cinsinden yükseklik
  units = "in"
)



















# 5. Yıllık Risk Resimleri oku ve ggplot nesnesine çevir
plot_list_yearly_risk <- lapply(SSP_LIST, function(s) {
  path <- here("outputs", s, "figures", "fig_yearly_risk.png")
  if (file.exists(path)) {
    img <- image_read(path)
    return(rasterGrob(img, interpolate = TRUE))
  }
  return(NULL)
})

# NULL olanları temizle
plot_list_yearly_risk <- plot_list_yearly_risk[!sapply(plot_list_yearly_risk, is.null)]

# 2. Patchwork ile birleştir (Örn: Yan yana 3 grafik)
combined_plot_yearly_risk <- wrap_plots(plot_list_yearly_risk) + 
  plot_layout(ncol = 2)

# 3. 300 DPI olarak kaydet
ggsave(
  filename = here(DIR_OUTPUT_CROSS, "combined_plot_yearly_risk.png"),
  plot = combined_plot_yearly_risk,
  dpi = 300,
  width = 15, # İnç cinsinden genişlik
  height = 10,  # İnç cinsinden yükseklik
  units = "in"
)
