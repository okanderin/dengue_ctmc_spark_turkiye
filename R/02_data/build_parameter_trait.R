source("R/01_setup/init.R")
source(here::here("R", "03_models", "parameter_functions.R"))

write_trait_param_template(file.path(DIR_PROCESSED, "trait_params.csv"))

