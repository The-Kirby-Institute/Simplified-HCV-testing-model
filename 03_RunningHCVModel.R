## run HCV model 

rm(list = ls()) 
# set memory limit 


library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
basePath <- getwd()
dataPath <- file.path(basePath, "01. Data/model input") 
Rcode <- file.path(basePath, "03. Code")
ResultsFolder <- file.path(basePath, "04. Output/Results")
currTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
project <- "Taiwanese MSM HCV model"
project_name <- "HCVModel"
projectFolder <- file.path(basePath)
source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))
source(file.path(Rcode, "/Functions/plotOptions.R"))


runSamples <- TRUE
saveAsBase <- FALSE  
runScenarios <- FALSE

# if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

# load project 
projectFile <- file.path(basePath,
                         paste0(project_name, ".rda"))
projectVars <- load(projectFile)
load(file.path(projectFolder, paste0(project_name, "cali",".rda")))

load(file.path(projectFolder, paste0(project_name, "param",".rda")))
load(file.path(projectFolder, paste0(project_name, "cost",".rda")))
# start clock 
tic <- proc.time()
# Run model on best estimates
bestResults <- HCVMSM(HCV,best_estimates, best_est_pop,
                        disease_progress,pop_array, dfList, fib,
                      end_Y =30, modelrun="UN", cost = cost, costflow = cost$flow)


# Run sampled parameter sets

HCV$numberSamples <- 1000
if (runSamples) {
  paramResults <- list()
  
  for (set in 1:HCV$numberSamples){
    
    
    paramResults[[set]] <- HCVMSM(HCV,
                                  Param_estimates[[set]], 
                                  Param_Pops[[set]], 
                                  Param_disease_progress[[set]], 
                                  Param_pop_array[[set]], 
                                  Param_dfList[[set]],
                                  Param_fib[[set]],
                                  end_Y = 30,modelrun="UN")
  }
 }
if (saveAsBase) {
  # Create testing folder and results file
  dir.create(file.path(ResultsFolder, "results_base"), 
             showWarnings = FALSE)
  resultsFile <- file.path(ResultsFolder, "results_base", 
                           "results_base.rda")
  
  save(list = c("bestResults", projectVars), 
       file = resultsFile) 
} else {
  dir.create(file.path(ResultsFolder, paste0("results_", currTime)), 
             showWarnings = FALSE)
  resultsFile <- file.path(ResultsFolder, paste0("results_", currTime), 
                           paste0("results_", currTime, ".rda"))
  
  save(list = c("bestResults", "paramResults",  projectVars), 
       file = resultsFile) 
}



#### Main scenarios & sensitivity #### 
# testing : change in the dfList 
# saving
# append all scenarios 
load(file.path(projectFolder, paste0(project_name, "Scenarios",".rda")))
tic <- proc.time()
scenarioResults <- list() 
for(i in 1: length(Scen_main)){
  
  scenarioResults[["Main"]][[names(Scen_main[i])]] <- 
    HCVMSM(HCV,best_estimates, best_est_pop, disease_progress,pop_array,
           Scen_main[[i]], fib, end_Y =30, modelrun="UN")
  
  scenarioResults[["tarpop"]][[names(Scen_main[i])]] <- 
    HCVMSM(HCV,best_estimates, best_est_pop, disease_progress,pop_array,
           Scen_targetpop[[i]], fib, end_Y =30, modelrun="UN")
  
  scenarioResults[["timeL"]][[names(Scen_main[i])]] <- 
    HCVMSM(HCV,best_estimates, best_est_pop, disease_progress,pop_array,
           Scen_timeL[[i]], fib, end_Y =30, modelrun="UN")
  
  scenarioResults[["timeU"]][[names(Scen_main[i])]] <- 
    HCVMSM(HCV,best_estimates, best_est_pop, disease_progress,pop_array,
           Scen_timeU[[i]], fib, end_Y =30, modelrun="UN")
  
  scenarioResults[["diagL"]][[names(Scen_main[i])]] <- 
    HCVMSM(HCV,best_estimates, best_est_pop, disease_progress,pop_array,
           Scen_diagL[[i]], fib, end_Y =30, modelrun="UN")
  
  scenarioResults[["diagU"]][[names(Scen_main[i])]] <- 
    HCVMSM(HCV,best_estimates, best_est_pop, disease_progress,pop_array,
           Scen_diagU[[i]], fib, end_Y =30, modelrun="UN")
  
  scenarioResults[["reinfL"]][[names(Scen_main[i])]] <- 
    HCVMSM(HCV,best_estimates, best_est_pop, disease_progress,pop_array,
           Scen_main[[i]], fib, end_Y =30, modelrun="UN", scenario = reinfL)
  
  scenarioResults[["reinfU"]][[names(Scen_main[i])]] <- 
    HCVMSM(HCV,best_estimates, best_est_pop, disease_progress,pop_array,
           Scen_main[[i]], fib, end_Y =30, modelrun="UN", scenario = reinfU)
  
}


save(scenarioResults,
     file = file.path(projectFolder,
                      paste0(project_name, "ScenResults", ".rda")))

rm(bestResults)
gc()
# uncertainty 
tic <- proc.time()
# main 
load(file.path(projectFolder, paste0(project_name, "Scenarios_mainpar",".rda")))

scenarioResultsPar_main <- list() 

for(i in c(SceName)){ 
  for(m in c(1: HCV$numberSamples)){
    scenarioResultsPar_main[[i]][[m]] <- 
      HCVMSM(HCV,
             Param_estimates[[m]], 
             Param_Pops[[m]], 
             Param_disease_progress[[m]], 
             Param_pop_array[[m]], 
             Scen_mainParam[[i]][[m]], Param_fib[[m]], end_Y =30, modelrun="UN")
    }
  }
save(scenarioResultsPar_main,
     file = file.path(projectFolder,
                      paste0(project_name, "ScenResultsPar_main", ".rda")))   

rm(scenarioResultsPar_main)
gc()

# targetpop
load(file.path(projectFolder, paste0(project_name, "Scenarios_tarpoppar",".rda")))

scenarioResultsPar_tarpop <- list()

for(i in c(SceName)){ 
  for(m in c(1: HCV$numberSamples)){    
    scenarioResultsPar_tarpop[[i]][[m]] <- 
      HCVMSM(HCV,
             Param_estimates[[m]], 
             Param_Pops[[m]], 
             Param_disease_progress[[m]], 
             Param_pop_array[[m]], 
             Scen_targetpopParam[[i]][[m]], Param_fib[[m]], end_Y =30, 
             modelrun="UN")
    }
}

save(scenarioResultsPar_tarpop,
     file = file.path(projectFolder,
                      paste0(project_name, "ScenResultsPar_tarpop", ".rda"))) 
rm(scenarioResultsPar_tarpop, Scen_targetpopParam)
gc()    

# time_lower
load(file.path(projectFolder, paste0(project_name, "Scenarios_timeLpar",".rda")))
scenarioResultsPar_timeL <- list()
for(i in c(SceName)){ 
  for(m in c(1: HCV$numberSamples)){      
    scenarioResultsPar_timeL[[i]][[m]] <- 
      HCVMSM(HCV,
             Param_estimates[[m]], 
             Param_Pops[[m]], 
             Param_disease_progress[[m]], 
             Param_pop_array[[m]], 
             Scen_timeLParam[[i]][[m]], Param_fib[[m]], end_Y =30, 
             modelrun="UN")
    }
}

save(scenarioResultsPar_timeL,
     file = file.path(projectFolder,
                      paste0(project_name, "ScenResultsPar_timeL", ".rda"))) 

rm(scenarioResultsPar_timeL, Scen_timeLParam)
gc()    
    
# timeU
load(file.path(projectFolder, paste0(project_name, "Scenarios_timeUpar",".rda")))

scenarioResultsPar_timeU <- list()
for(i in c(SceName)){ 
  for(m in c(1: HCV$numberSamples)){ 
    scenarioResultsPar_timeU[[i]][[m]] <- 
      HCVMSM(HCV,
             Param_estimates[[m]], 
             Param_Pops[[m]], 
             Param_disease_progress[[m]], 
             Param_pop_array[[m]], 
             Scen_timeUParam[[i]][[m]], Param_fib[[m]], end_Y =30, 
             modelrun="UN") 
    }
  }
save(scenarioResultsPar_timeU,
     file = file.path(projectFolder,
                      paste0(project_name, "ScenResultsPar_timeU", ".rda"))) 

rm(scenarioResultsPar_timeU, Scen_timeUParam)
gc()

# diagL
load(file.path(projectFolder, paste0(project_name, "Scenarios_diagLpar",".rda")))
scenarioResultsPar_diagL <- list()
for(i in c(SceName)){ 
  for(m in c(1: HCV$numberSamples)){     
    scenarioResultsPar_diagL[[i]][[m]] <- 
      HCVMSM(HCV,
             Param_estimates[[m]], 
             Param_Pops[[m]], 
             Param_disease_progress[[m]], 
             Param_pop_array[[m]], 
             Scen_diagLParam[[i]][[m]], Param_fib[[m]], end_Y =30, 
             modelrun="UN")
    }
  }
save(scenarioResultsPar_diagL,
     file = file.path(projectFolder,
                      paste0(project_name, "ScenResultsPar_diagL", ".rda"))) 

rm(scenarioResultsPar_diagL, Scen_diagLParam)
gc()

# diagU
load(file.path(projectFolder, paste0(project_name, "Scenarios_diagUpar",".rda")))
scenarioResultsPar_diagU <- list()
for(i in c(SceName)){ 
  for(m in c(1: HCV$numberSamples)){    

    scenarioResultsPar_diagU[[i]][[m]] <- 
      HCVMSM(HCV,
             Param_estimates[[m]], 
             Param_Pops[[m]], 
             Param_disease_progress[[m]], 
             Param_pop_array[[m]], 
             Scen_diagUParam[[i]][[m]], Param_fib[[m]], end_Y =30, 
             modelrun="UN")
    }
  }
save( scenarioResultsPar_diagU,
     file = file.path(projectFolder,
                      paste0(project_name, "ScenResultsPar_diagU", ".rda"))) 

rm( scenarioResultsPar_diagU, Scen_diagUParam)
gc()  

# reinfL
load(file.path(projectFolder, paste0(project_name, "Scenarios",".rda")))
scenarioResultsPar_reinfL <- list()
for(i in c(SceName)){ 
  for(m in c(1: HCV$numberSamples)){  
    
    scenarioResultsPar_reinfL[[i]][[m]] <- 
      HCVMSM(HCV,
             Param_estimates[[m]], 
             Param_Pops[[m]], 
             Param_disease_progress[[m]], 
             Param_pop_array[[m]], 
             Scen_mainParam[[i]][[m]], Param_fib[[m]], end_Y =30, modelrun="UN",
             scenario = reinfL)
    }
}

save(scenarioResultsPar_reinfL,
      file = file.path(projectFolder,
                       paste0(project_name, "ScenResultsPar_reinfL", ".rda"))) 


gc() 
    
#reinfU
scenarioResultsPar_reinfU <- list()
for(i in c(SceName)){ 
  for(m in c(1: HCV$numberSamples)){ 
    scenarioResultsPar_reinfU[[i]][[m]] <- 
      HCVMSM(HCV,
             Param_estimates[[m]], 
             Param_Pops[[m]], 
             Param_disease_progress[[m]], 
             Param_pop_array[[m]], 
             Scen_mainParam[[i]][[m]], Param_fib[[m]], end_Y =30, modelrun="UN",
             scenario = reinfU)
    
    }
  }
save(scenarioResultsPar_reinfU,
     file = file.path(projectFolder,
                      paste0(project_name, "ScenResultsPar_reinfU", ".rda"))) 


gc()


toc <- proc.time() - tic 
toc



