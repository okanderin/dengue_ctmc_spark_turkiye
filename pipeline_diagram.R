# =========================================================
# Publication-quality project structure diagram for r_project_tez
# Based on TREE /F output
# =========================================================

# ---------- Packages ----------
req_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ---------- Graph ----------
gr <- grViz("
digraph r_project_tez_structure {

  graph [
    layout = dot,
    rankdir = TB,
    splines = ortho,
    nodesep = 0.35,
    ranksep = 0.55,
    pad = 0.25,
    bgcolor = 'white',
    fontname = 'Helvetica'
  ]

  node [
    shape = box,
    style = 'rounded,filled',
    fontname = 'Helvetica',
    fontsize = 11,
    margin = '0.12,0.08',
    color = '#5B6470',
    penwidth = 1.0,
    fillcolor = '#F8FAFC',
    fontcolor = '#1F2937'
  ]

  edge [
    color = '#6B7280',
    penwidth = 1.1,
    arrowsize = 0.75
  ]

  # ======================================================
  # ROOT
  # ======================================================
  root [
    label = 'r_project_tez\\nProject root',
    shape = box,
    style = 'rounded,filled,bold',
    fillcolor = '#111827',
    color = '#111827',
    fontcolor = 'white',
    fontsize = 16,
    margin = '0.18,0.12'
  ]

  readme [label = 'README.md']
  rproj  [label = 'r_project_tez.Rproj']
  exfile [label = 'örnek_for_sspler_icin.R']

  root -> readme
  root -> rproj
  root -> exfile

  # ======================================================
  # DATA RAW
  # ======================================================
  subgraph cluster_raw {
    label = 'Raw inputs';
    labelloc = 't';
    fontsize = 14;
    fontname = 'Helvetica';
    color = '#93C5FD';
    penwidth = 1.3;
    style = 'rounded';
    bgcolor = '#EFF6FF';

    data_raw [label = 'data_raw', fillcolor = '#DBEAFE', color = '#60A5FA', style='rounded,filled,bold']

    raw_climate [label = 'climate\\nCMIP6 + historical NetCDF']
    raw_geo     [label = 'geography\\ndistrict shapefiles']
    raw_meta    [label = 'metadata\\nimportation metadata']
    raw_pop     [label = 'population\\nGBD + TÜİK + tourism']
    raw_vector  [label = 'vector\\npresence + literature parameters']

    cmip6      [label = 'cmip6\\nssp126 / ssp245 / ssp585']
    hist       [label = 'cmip6_hist_1981_2014']
    geo_shp    [label = 'districts_shapefile\\nGADM layers']
    pop_gbd    [label = 'GBD_2023_DATA']
    pop_proj   [label = 'projections']
    pop_tour   [label = 'turizm_verileri']
    pop_turk   [label = 'turkstat']

    data_raw -> raw_climate
    data_raw -> raw_geo
    data_raw -> raw_meta
    data_raw -> raw_pop
    data_raw -> raw_vector

    raw_climate -> cmip6
    raw_climate -> hist
    raw_geo -> geo_shp
    raw_pop -> pop_gbd
    raw_pop -> pop_proj
    raw_pop -> pop_tour
    raw_pop -> pop_turk
  }

  # ======================================================
  # DATA INTERIM
  # ======================================================
  subgraph cluster_interim {
    label = 'Interim climate objects';
    labelloc = 't';
    fontsize = 14;
    color = '#86EFAC';
    penwidth = 1.3;
    style = 'rounded';
    bgcolor = '#F0FDF4';

    data_interim [label = 'data_interim', fillcolor = '#DCFCE7', color = '#4ADE80', style='rounded,filled,bold']

    interim126 [label = 'ssp126\\nclimate_district_monthly_2015_2100.rds']
    interim245 [label = 'ssp245\\nclimate_district_monthly_2015_2100.rds']
    interim585 [label = 'ssp585\\nclimate_district_monthly_2015_2100.rds']

    data_interim -> interim126
    data_interim -> interim245
    data_interim -> interim585
  }

  # ======================================================
  # DATA PROCESSED
  # ======================================================
  subgraph cluster_processed {
    label = 'Processed analytical datasets';
    labelloc = 't';
    fontsize = 14;
    color = '#C4B5FD';
    penwidth = 1.3;
    style = 'rounded';
    bgcolor = '#F5F3FF';

    data_processed [label = 'data_processed', fillcolor = '#EDE9FE', color = '#A78BFA', style='rounded,filled,bold']

    proc_shared [label = 'Shared processed files\\nbias correction, population, seasonality,\\nspecies, traits, travel weights']
    proc126     [label = 'ssp126\\nclimate + importation products']
    proc245     [label = 'ssp245\\nclimate + importation products']
    proc585     [label = 'ssp585\\nclimate + importation products']

    proc126_2   [label = '2024_2100\\nmonthly/yearly importation']
    proc245_2   [label = '2024_2100\\nmonthly/yearly importation']
    proc585_2   [label = '2024_2100\\nmonthly/yearly importation']

    data_processed -> proc_shared
    data_processed -> proc126
    data_processed -> proc245
    data_processed -> proc585

    proc126 -> proc126_2
    proc245 -> proc245_2
    proc585 -> proc585_2
  }

  # ======================================================
  # LOGS
  # ======================================================
  subgraph cluster_logs {
    label = 'Execution logs';
    labelloc = 't';
    fontsize = 14;
    color = '#F9A8D4';
    penwidth = 1.3;
    style = 'rounded';
    bgcolor = '#FDF2F8';

    logs [label = 'logs', fillcolor = '#FCE7F3', color = '#F472B6', style='rounded,filled,bold']
    log_types [label = 'Pipeline logs\\nimport_raw / clean_interim /\\nbuild_population / travel_weights /\\nimportation_pressure / SSP-specific logs']
    logs -> log_types
  }

  # ======================================================
  # OUTPUTS
  # ======================================================
  subgraph cluster_outputs {
    label = 'Model outputs';
    labelloc = 't';
    fontsize = 14;
    color = '#FCD34D';
    penwidth = 1.3;
    style = 'rounded';
    bgcolor = '#FFFBEB';

    outputs [label = 'outputs', fillcolor = '#FEF3C7', color = '#F59E0B', style='rounded,filled,bold']

    cross_scenario [label = 'cross_scenario']
    out126 [label = 'ssp126']
    out245 [label = 'ssp245']
    out585 [label = 'ssp585']

    outputs -> cross_scenario
    outputs -> out126
    outputs -> out245
    outputs -> out585

    # ----- SSP126 detail block -----
    out126_diag  [label = 'diagnostics\\nMC validation, R0 checks']
    out126_fig   [label = 'figures\\nheatmaps, seasonality,\\nvalidation, yearly risk']
    out126_mod   [label = 'model_results\\nmonthly / yearly / horizon CTMC']
    out126_rain  [label = 'rainfall_CSI\\nCSI and seasonal modulation']
    out126_sens  [label = 'sensitivity\\nctmc_mc + importation']
    out126_sim   [label = 'simulation\\nrep1000 outputs']
    out126_tab   [label = 'tables\\nsummary and regression tables']

    out126 -> out126_diag
    out126 -> out126_fig
    out126 -> out126_mod
    out126 -> out126_rain
    out126 -> out126_sens
    out126 -> out126_sim
    out126 -> out126_tab

    # ----- SSP245 and SSP585 summarized -----
    out245_sum [label = 'parallel SSP245 structure\\n(diagnostics / figures / model_results /\\nsensitivity / simulation / tables)']
    out585_sum [label = 'parallel SSP585 structure\\n(diagnostics / figures / model_results /\\nsensitivity / simulation / tables)']

    out245 -> out245_sum
    out585 -> out585_sum
  }

  # ======================================================
  # FLOW RELATIONS
  # ======================================================
  root -> data_raw
  root -> data_interim
  root -> data_processed
  root -> logs
  root -> outputs

  raw_climate -> data_interim [ltail=cluster_raw, lhead=cluster_interim]
  data_interim -> data_processed [ltail=cluster_interim, lhead=cluster_processed]
  data_processed -> outputs [ltail=cluster_processed, lhead=cluster_outputs]
  data_processed -> logs
}
")

# Show in viewer
gr

# ---------- Export ----------
svg_txt <- export_svg(gr)
writeLines(svg_txt, 'r_project_tez_structure_publication.svg')
rsvg_png(charToRaw(svg_txt), file = 'r_project_tez_structure_publication.png', width = 7200, height = 9600)

# Optional PDF export (requires cairo-capable setup on some systems)
# rsvg_pdf(charToRaw(svg_txt), file = 'r_project_tez_structure_publication.pdf')