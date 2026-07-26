# ==========================================================
# Kavramsal Çerçeve Diyagramı (Sadece Stokastik Spark Fazı)
# - DiagrammeR ile çizim
# - DiagrammeRsvg + rsvg ile SVG/PNG çıktı
# - Akış: Girdiler -> P_est + q_import_month -> p_month_major
#         -> p_year -> p_horizon (2025–2075)
# ==========================================================

# ---- Paketler ----
req_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ---- Dosya adları (çıktılar) ----
out_svg <- "kavramsal_cerceve_spark_only.svg"
out_png <- "kavramsal_cerceve_spark_only.png"

# ---- Diyagram ----
kavram_spark_only <- grViz("
digraph spark_only {
  graph [layout = dot, rankdir = LR, nodesep = 0.60, ranksep = 0.95]

  node [shape=rectangle, style=\"filled,rounded\", fontname=\"Times-Roman\", fontsize=12, penwidth=2.5]
  edge [fontname=\"Times-Roman\", fontsize=10, color=\"#404040\", penwidth=1.5, fontcolor=\"#303030\", arrowsize=0.8]

  // ====== GİRDİ HAVUZU ======
  subgraph cluster_inputs {
    label = \"Girdi ve Parametre Havuzu\";
    style = dashed; color = gray55; margin = 18;

    import [label = \"İthal Baskı\\nq_import, λ_import\\n(Yolcu + Kuresel insidans)\", fillcolor = \"#D1E8E2\"]
    climate [label = \"Iklim & Demografi\\nAylık T, RH (CMIP6 SSP2-4.5)\\nN_h (TUIK/IIASA)\", fillcolor = \"#D1E8E2\"]
    traits  [label = \"Biyolojik/Traits\\na(T), EIP(T), μ_v(T,RH), β, m\\n(Spp: aegypti/albopictus)\", fillcolor = \"#FFF2CC\"]

    {rank=same; import; climate; traits}
  }

  // ====== ÇEKİRDEK CTMC / HESAP ======
  subgraph cluster_core {
    label = \"Stokastik (CTMC)\"; style = rounded; color = \"#8EC07C\"; penwidth=2;

    pest [label = \"Yerel Tutunma\\nP_est = 1 - q_extinction\\n(BD: λ_local, μ = γ·I)\", fillcolor = \"#D9EAD3\"]
    qimp [label = \"O Ay ≥1 İthal Giriş\\nq_import,month = 1 - e^{-Λ_import}\", fillcolor = \"#E6F5F3\"]
    pmonth [label = \"Aylık Major Olasılık\\np_month,major = q_import,month × P_est\", fillcolor = \"#CCE3F5\", penwidth=3]
  }

  // ====== ZAMAN OLÇEKLERİ / BİRLESTİRME ======
  yearly   [label = \"Yıllık Risk\\np_year = 1 - Π_m (1 - p_month,major)\", fillcolor = \"#E1D5E7\"]
  horizon  [label = \"Ufuk (2025–2075)\\np_horizon = 1 - Π_t (1 - p_month,major)\", fillcolor = \"#E1D5E7\", penwidth=3]

  // ====== OKLAR ======
  import  -> qimp   [label = \"Λ_import, q_import,month\"]
  traits  -> pest   [label = \"λ_local(T,RH), γ\" color=\"#E67E22\" fontcolor=\"#E67E22\"]
  climate -> pest   [style=dotted, label = \"T, RH\"]
  climate -> traits [style=dotted, label = \"Parametrelemede kullanılır\"]
  traits  -> qimp   [style=dotted, label = \"(opsiyonel) M_climate(t)\" color=\"#888888\" fontcolor=\"#666666\"]

  pest -> pmonth
  qimp -> pmonth

  pmonth -> yearly [label = \"Aylık -> Yıllık Bernoulli\"]
  yearly -> horizon [label = \"Yıllar Boyunca Birikim (2025–2075)\"]
}
")

# ---- SVG/PNG kaydet ----
svg_txt <- export_svg(kavram_spark_only)
# PNG genişlik/yükseklik: poster/tez için yüksek çözünürlük
rsvg_svg(charToRaw(svg_txt), file = out_svg)
rsvg_png(charToRaw(svg_txt), file = out_png, width = 3600, height = 900)

# Ekrana da çiz
kavram_spark_only



# ==========================================================
# Kavramsal Çerçeve Diyagramı (Sadece Stokastik Spark Fazı)
# - Yıllık risk ve ufuk riski kutuları aylık riskin altına alındı
# - Manuel yerleşim: neato + pos
# ==========================================================

# ---- Paketler ----
req_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ---- Dosya adları (çıktılar) ----
out_svg <- "kavramsal_cerceve_spark_only_verticalrisk.svg"
out_png <- "kavramsal_cerceve_spark_only_verticalrisk.png"

# ---- Diyagram ----
kavram_spark_only <- grViz("
digraph spark_only {

  graph [
    layout = neato,
    overlap = false,
    splines = ortho,
    outputorder = edgesfirst,
    pad = 0.4
  ]

  node [
    shape = rectangle,
    style = \"filled,rounded\",
    fontname = \"Times-Roman\",
    fontsize = 12,
    penwidth = 2.5,
    width = 2.4,
    height = 0.9
  ]

  edge [
    fontname = \"Times-Roman\",
    fontsize = 10,
    color = \"#404040\",
    penwidth = 1.5,
    fontcolor = \"#303030\",
    arrowsize = 0.8
  ]

  // ====== GİRDİ HAVUZU ======
  subgraph cluster_inputs {
    label = \"Girdi ve Parametre Havuzu\"
    style = dashed
    color = gray55
    margin = 20

    import [
      label = \"İthal Baskı\\nq_import, λ_import\\n(Yolcu + Küresel insidans)\",
      fillcolor = \"#D1E8E2\",
      pos = \"0,3!\"
    ]

    climate [
      label = \"İklim & Demografi\\nAylık T, RH (CMIP6 SSP2-4.5)\\nN_h (TÜİK/IIASA)\",
      fillcolor = \"#D1E8E2\",
      pos = \"0,-2!\"
    ]

    traits [
      label = \"Biyolojik/Traits\\na(T), EIP(T), μ_v(T,RH), β, m\\n(Spp: aegypti/albopictus)\",
      fillcolor = \"#FFF2CC\",
      pos = \"0,0.5!\"
    ]
  }

  // ====== ÇEKİRDEK CTMC / HESAP ======
  subgraph cluster_core {
    label = \"Stokastik (CTMC)\"
    style = rounded
    color = \"#8EC07C\"
    penwidth = 2
    margin = 20

    qimp [
      label = \"O Ay ≥1 İthal Giriş\\nq_import,month = 1 - e^{-Λ_import}\",
      fillcolor = \"#E6F5F3\",
      pos = \"5,2.2!\"
    ]

    pest [
      label = \"Yerel Tutunma\\nP_est = 1 - q_extinction\\n(BD: λ_local, μ = γ·I)\",
      fillcolor = \"#D9EAD3\",
      pos = \"5,0!\"
    ]

    pmonth [
      label = \"Aylık Major Olasılık\\np_month,major = q_import,month × P_est\",
      fillcolor = \"#CCE3F5\",
      penwidth = 3,
      pos = \"9,1.1!\"
    ]
  }

  // ====== ZAMAN ÖLÇEKLERİ / BİRLEŞTİRME ======
  yearly [
    label = \"Yıllık Risk\\np_year = 1 - Π_m (1 - p_month,major)\",
    fillcolor = \"#E1D5E7\",
    pos = \"9,-1.3!\"
  ]

  horizon [
    label = \"Ufuk (2025–2075)\\np_horizon = 1 - Π_t (1 - p_month,major)\",
    fillcolor = \"#E1D5E7\",
    penwidth = 3,
    pos = \"9,-3.7!\"
  ]

  // ====== OKLAR ======
  import  -> qimp   [label = \"Λ_import, q_import,month\"]
  traits  -> pest   [label = \"λ_local(T,RH), γ\", color = \"#E67E22\", fontcolor = \"#E67E22\"]
  climate -> pest   [style = dotted, label = \"T, RH\"]
  climate -> traits [style = dotted, label = \"Parametrelemede kullanılır\"]
  traits  -> qimp   [style = dotted, label = \"(opsiyonel) M_climate(t)\", color = \"#888888\", fontcolor = \"#666666\"]

  pest   -> pmonth
  qimp   -> pmonth
  pmonth -> yearly  [label = \"Aylık -> Yıllık Birikim\"]
  yearly -> horizon [label = \"Yıllar Boyunca Birikim (2025–2075)\"]
}
")

# ---- SVG/PNG kaydet ----
svg_txt <- export_svg(kavram_spark_only)

rsvg_svg(charToRaw(svg_txt), file = out_svg)
rsvg_png(charToRaw(svg_txt), file = out_png, width = 3200, height = 1800)

# ---- Ekrana çiz ----
kavram_spark_only




