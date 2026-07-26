# ==========================================================
# Kavramsal Çerçeve Diyagramı — DİKEY/KARE düzen (Word okunurluğu)
# - rankdir=TB: akış yukarıdan aşağı
# - Girdiler üstte yan yana, CTMC ortada, sonuç zinciri altta DİKEY tek kolon
# - Oran ~1.5 (3.6 yerine) -> Word'de kutular büyük ve okunur
# - SVG'yi Word 365'e gömün (Ekle > Resimler > Bu Cihaz)
# ==========================================================

req_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

out_svg <- "kavramsal_cerceve_kare.svg"
out_png <- "kavramsal_cerceve_kare.png"

kavram <- grViz("
digraph spark {
  graph [layout = dot, rankdir = TB, nodesep = 0.45, ranksep = 0.70, newrank = true]
  node [shape=rectangle, style=\"filled,rounded\", fontname=\"Times-Roman\", fontsize=14, penwidth=2.5, margin=\"0.20,0.12\"]
  edge [fontname=\"Times-Roman\", fontsize=12, color=\"#404040\", penwidth=1.6, fontcolor=\"#303030\", arrowsize=0.85]

  // ÜST: Girdiler yan yana
  subgraph cluster_inputs {
    label = \"Girdi ve Parametre Havuzu\"; labelloc=t;
    style = dashed; color = gray55; margin = 14; fontname=\"Times-Roman\"; fontsize=14;
    import  [label = \"İthal Baskı\\nq_import, λ_import\\n(Yolcu + Küresel insidans)\", fillcolor = \"#D1E8E2\"]
    traits  [label = \"Biyolojik / Traits\\na(T), EIP(T), μ_v(T,RH), β, m\\n(aegypti / albopictus)\", fillcolor = \"#FFF2CC\"]
    climate [label = \"İklim & Demografi\\nAylık T, RH (CMIP6 SSP2-4.5)\\nN_h (TÜİK/IIASA)\", fillcolor = \"#D1E8E2\"]
    {rank=same; import; traits; climate}
  }

  // ORTA: CTMC ikili
  subgraph cluster_core {
    label = \"Stokastik (CTMC)\"; labelloc=t;
    style = rounded; color = \"#8EC07C\"; penwidth=2; margin=14; fontname=\"Times-Roman\"; fontsize=14;
    qimp [label = \"O Ay ≥1 İthal Giriş\\nq_import,month = 1 − e^(−Λ_import)\", fillcolor = \"#E6F5F3\"]
    pest [label = \"Yerel Tutunma\\nP_est = 1 − q_extinction\\n(BD: λ_local, μ = γ·I)\", fillcolor = \"#D9EAD3\"]
    {rank=same; qimp; pest}
  }

  // ALT: sonuç zinciri DİKEY tek kolon
  pmonth  [label = \"Aylık Major Olasılık\\np_month,major = q_import,month × P_est\", fillcolor = \"#CCE3F5\", penwidth=3]
  yearly  [label = \"Yıllık Risk\\np_year = 1 − Π_m (1 − p_month,major)\", fillcolor = \"#E1D5E7\"]
  horizon [label = \"Ufuk (2025–2075)\\np_horizon = 1 − Π_t (1 − p_month,major)\", fillcolor = \"#E1D5E7\", penwidth=3]

  // OKLAR
  climate -> traits [style=dotted, label = \"Parametreleme\"]
  climate -> pest   [style=dotted, label = \"T, RH\"]
  import  -> qimp   [label = \"Λ_import\"]
  traits  -> pest   [label = \"λ_local, γ\", color=\"#E67E22\", fontcolor=\"#E67E22\"]
  traits  -> qimp   [style=dotted, label = \"(ops.) M_climate(t)\", color=\"#888888\", fontcolor=\"#666666\"]

  qimp -> pmonth
  pest -> pmonth
  pmonth -> yearly  [label = \"Aylık → Yıllık birikim\"]
  yearly -> horizon [label = \"2025–2075 birikim\"]
}
")

# ---- Kaydet ----
svg_txt <- export_svg(kavram)
rsvg_svg(charToRaw(svg_txt), file = out_svg)   # <- Word'e bunu gömün
rsvg_png(charToRaw(svg_txt), file = out_png, width = 2400, height = 1568)  # yedek yüksek çöz. PNG

kavram