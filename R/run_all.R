
------------------------
  # R/run_all.R
-----------------------




source("R/01_setup/init.R")
source("R/03_models/parameter_functions.R")

source("R/04_results/make_delta_metrics.R")
source("R/04_results/derive_core_metrics.R")
source("R/04_results/core_metrics_calc.R")

--------------------------
#duyarlılık analizi
--------------------------
source("R/04_results/run_m_sensitivity.R")

source("R/04_results/make_plots.R")