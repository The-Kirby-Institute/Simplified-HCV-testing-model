#### run scripts #### 

# after running {scenario set up_20260122_scaledscaleup.R}
rm(list = ls())
gc()
source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/02_3. Uncertainty.R", echo = TRUE)
rm(list = ls())
gc()
source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/02_3_1. Scenario_uncertainty.R", echo = TRUE)
rm(list = ls())
gc()
source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/02_4. Param_sim.R", echo = TRUE)
rm(list = ls())
gc()
source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/02_4. Param_sim_scenario.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/Resultsummary.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/Resultsummary_sensitivity.R", echo = TRUE)
rm(list = ls())
gc()

source("~/Projects/Simplified-HCV-testing-model/Projects/POC_AU/Res_aggregate_sensitivity.R", echo = TRUE)
rm(list = ls())
gc()

