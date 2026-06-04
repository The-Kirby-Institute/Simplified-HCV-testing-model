#### run scripts #### 

# after calibrated to 2021
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/calibration_pipeline.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/scenarios_simulation_testing number.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/02_3_1. Scenario_uncertainty.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/02_4. Param_sim_scenario.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/Resultsummary.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/Res_aggregate_sensitivity.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/final_figs.R", echo = TRUE)