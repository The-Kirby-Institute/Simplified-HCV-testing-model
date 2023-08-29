# simulating scenarios under uptaking PrEP and HIVD 
rm(list = ls())

# check work directory
library(here)
here()

# load libraries  
library(dplyr)
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library(doParallel)

# we specify the number of cores/workers we want to use
registerDoParallel(cores = detectCores() - 1)

#file path of "TWHCV-model" project
codepath <- file.path(here() %>% dirname(), 'TWHCV-model/03. Code/Functions')

DataFolder <- file.path(here(), "01. DATA/model input")

# Rda file path 
# load the .rda file of base estimate 
rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")



projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile)

load(file.path(rdapath , paste0("HCVModel", "cali",".rda")))

load(file.path(here("02. Output/RDA"), "cost.rda"))

load(file.path(here("02. Output/RDA"), "Param_dfExt.rda"))

load(file.path(here("02. Output/RDA"), "mainS_PrEPHIV.rda"))

load(file.path(here("02. Output/RDA"), "mainansis.rda"))

load(file.path(here("02. Output/RDA"), "toft.rda"))

# output file path 
# dir.create("02. Output") # create subdircetory 
# dir.create("02. Output/RDA") 
# dir.create("02. Output/Results")
ResultsFolder <- file.path(here(), "02. Output/Results")
outputdt <- here("02. Output/RDA")

# source 
source(file.path(codepath, "HCV_model.R"))

source(file.path(codepath, "plotFunctions.R"))

source(file.path(codepath, "plotOptions.R"))

source(file.path(codepath, "Scenarios_wrangling.R"))

currTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

runSamples <- TRUE

saveAsBase <- FALSE  

runScenarios <- FALSE 

simY <- 100

# Run model on best estimates
# costTrim for scenarions 
CostScen <- list()
CostScen[["POC_antibody"]] <- costTrim$POCab
CostScen[["DBS"]] <- costTrim$POCab
CostScen[["Reflex_RNA"]] <- costTrim$Reflex
CostScen[["POC_RNA"]] <- costTrim$POCRNA


CostScen_par <- list()

for(set in 1: HCV$numberSamples){ 
  CostScen_par[[set]] <- CostScen
  
  }

for(set in 1:HCV$numberSamples){ 
  CostScen_par[[set]][["POC_antibody"]] <- Param_costList[[set]]$POCab
  CostScen_par[[set]][["DBS"]] <- Param_costList[[set]]$POCab
  CostScen_par[[set]][["Reflex_RNA"]] <- Param_costList[[set]]$Reflex
  CostScen_par[[set]][["POC_RNA"]] <- Param_costList[[set]]$POCRNA
  }



ScenResults <- list()

for(i in names(CEA_main)){ 
  
  ScenResults[[i]] <- HCVMSM(HCV,best_estimatesOff, best_est_pop,
                        disease_progress,pop_array, CEA_main[[i]], fib,
                        end_Y = simY, modelrun="UN", cost = costTrim, 
                        costflow = CostScen[[i]])
    }

# for scenarios related to uptake PrEP and HIV diagnosis rate
ScenResults_PrEPHIV <- list()

for(i in names(pop_arrayPrEPHIV)){
  
  for(m in names(CEA_main)){
    
    ScenResults_PrEPHIV[[i]][[m]] <- HCVMSM(HCV,best_estimatesOff, best_est_pop,
                                       disease_progress,pop_arrayPrEPHIV[[i]], 
                                       CEA_main[[m]], fib,
                                       end_Y = simY, modelrun="UN", cost = costTrim, 
                                       costflow = CostScen[[m]])
    }
  }

# Run sampled parameter sets
HCV$numberSamples <- 1000

load(file.path(rdapath , paste0("HCVModel", "param",".rda")))


# trim the time points 
trimpt <- 1000

for(set in 1:HCV$numberSamples){ 
  
  Param_pop_array[[set]] <- Param_pop_array[[set]][ , ,c(1:trimpt)]
  
}
### 

tic <- proc.time()

if (runSamples) {
  
  ScenparamResults <- list()

    for(m in names(CEA_main)){
      for (set in 1:HCV$numberSamples){
        
        ScenparamResults[[m]][[set]] <- HCVMSM(HCV,
                                                    Param_estimatesOff[[set]], 
                                                    Param_Pops[[set]], 
                                                    Param_disease_progress[[set]], 
                                                    Param_pop_array[[set]], 
                                                    CEA_mainParam[[m]][[set]],
                                                    Param_fib[[set]],
                                                    end_Y = simY,modelrun="UN",
                                                    cost = Param_costList[[set]], 
                                                    costflow = CostScen_par[[set]][[m]])
        }
      }
    }


toc <- proc.time() - tic

toc

save(ScenResults, ScenparamResults,
     file = file.path(outputdt, "simScen.rda"))

save(ScenResults_PrEPHIV,
     file = file.path(outputdt, "simScen_PrEPHIV.rda"))

rm(ScenResults, ScenparamResults, ScenResults_PrEPHIV)


# parameter sets for scenarios with uptake PrEP and HIVD 
ScenparamResults_PrEP <- list()
  
ScenparamResults_HIVD <- list()
  
ScenparamResults_PrEPnHIVD <- list()

tic <- proc.time()

for(m in names(CEA_main)){
  for (set in 1:HCV$numberSamples){
      
      ScenparamResults_PrEP[[m]][[set]] <- HCVMSM(HCV,
                                       Param_estimatesOff[[set]], 
                                       Param_Pops[[set]], 
                                       Param_disease_progress[[set]], 
                                       Param_pop_arrayPrEPHIV[["PrEP"]][[set]], 
                                       CEA_mainParam[[m]][[set]],
                                       Param_fib[[set]],
                                       end_Y = simY,modelrun="UN",
                                       cost = Param_costList[[set]], 
                                       costflow = CostScen_par[[set]][[m]])
      }
}

toc <- proc.time() - tic 

toc 

save(ScenparamResults_PrEP,
     file = file.path(outputdt, "simScen_param_PrEP.rda"))

rm(ScenparamResults_PrEP)

tic <- proc.time()

for(m in names(CEA_main)){
  for (set in 1:HCV$numberSamples){
    
    ScenparamResults_HIVD[[m]][[set]] <- HCVMSM(HCV,
                                                Param_estimatesOff[[set]], 
                                                Param_Pops[[set]], 
                                                Param_disease_progress[[set]], 
                                                Param_pop_arrayPrEPHIV[["HIVD"]][[set]], 
                                                CEA_mainParam[[m]][[set]],
                                                Param_fib[[set]],
                                                end_Y = simY,modelrun="UN",
                                                cost = Param_costList[[set]], 
                                                costflow = CostScen_par[[set]][[m]])
  }
}

toc <- proc.time() - tic 

toc 

save(ScenparamResults_HIVD,
     file = file.path(outputdt, "simScen_param_HIVD.rda"))

rm(ScenparamResults_HIVD)

tic <- proc.time()

for(m in names(CEA_main)){
  for (set in 1:HCV$numberSamples){
    
    ScenparamResults_PrEPnHIVD[[m]][[set]] <- HCVMSM(HCV,
                                                Param_estimatesOff[[set]], 
                                                Param_Pops[[set]], 
                                                Param_disease_progress[[set]], 
                                                Param_pop_arrayPrEPHIV[["PrEPnHIVD"]][[set]], 
                                                CEA_mainParam[[m]][[set]],
                                                Param_fib[[set]],
                                                end_Y = simY,modelrun="UN",
                                                cost = Param_costList[[set]], 
                                                costflow = CostScen_par[[set]][[m]])
  }
}

toc <- proc.time() - tic 

toc 

save(ScenparamResults_PrEPnHIVD,
     file = file.path(outputdt, "simScen_param_PrEPnHIVD.rda"))

rm(ScenparamResults_PrEPnHIVD)



