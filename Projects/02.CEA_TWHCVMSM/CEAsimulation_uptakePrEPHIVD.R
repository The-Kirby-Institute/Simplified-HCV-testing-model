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

codep <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/TWHCV-model"
com_codeP <- "/Users/jjwu/Documents/Simplified-HCV-testing-model/03. Code"
epidatapath <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Taiwan-MSM-HCV-model"
#file path of "TWHCV-model" project
Rcode <- file.path(codep, '03. Code')
scenariopath <- file.path(epidatapath, '01. DATA/model input/Scenarios')

# save and cost data folder 
dtp <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/02.CEA_TWHCVMSM"
DataFolder <- file.path(dtp, "01. DATA/model input")
outputdt <- file.path(dtp, "02. Output/RDA")


projectFile <- file.path(epidatapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile)

load(file.path(epidatapath , paste0("HCVModel", "cali",".rda")))

load(file.path(outputdt, "cost.rda"))

load(file.path(outputdt, "Param_dfExt.rda"))

load(file.path(outputdt, "mainS_PrEPHIV.rda"))

load(file.path(outputdt, "mainansis.rda"))

load(file.path(outputdt, "toft.rda"))

# output file path 
# dir.create("02. Output") # create subdircetory 
# dir.create("02. Output/RDA") 
# dir.create("02. Output/Results")
ResultsFolder <- file.path(outputdt, "02. Output/Results")


# source 
source(file.path(Rcode, "Functions/HCV_model.R"))

source(file.path(com_codeP, "Functions/plotFunctions.R"))

source(file.path(com_codeP, "Functions/plotOptions.R"))

source(file.path(com_codeP, "Functions/Scenarios_wrangling.R"))

currTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

runSamples <- TRUE

saveAsBase <- FALSE  

runScenarios <- FALSE 

simY <- 100

# Run model on best estimates
# costTrim for scenarions 
CostScen <- list()
CostScen[["POC_antibody"]] <- costTrim$POCab
CostScen[["Reflex_RNA"]] <- costTrim$Reflex
CostScen[["POC_RNA"]] <- costTrim$POCRNA


CostScen_par <- list()

for(set in 1: HCV$numberSamples){ 
  CostScen_par[[set]] <- CostScen
  
  }

for(set in 1:HCV$numberSamples){ 
  CostScen_par[[set]][["POC_antibody"]] <- Param_costList[[set]]$POCab
  CostScen_par[[set]][["Reflex_RNA"]] <- Param_costList[[set]]$Reflex
  CostScen_par[[set]][["POC_RNA"]] <- Param_costList[[set]]$POCRNA
  }

# remove "DBS" list

CEA_main <- CEA_main[names(CEA_main)!= "DBS"]

CEA_mainParam <- CEA_mainParam[names(CEA_mainParam)!= "DBS"]
# uptake HCV testing  
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

load(file.path(epidatapath , paste0("HCVModel", "param",".rda")))


# trim the time points 
trimpt <- 1000

for(set in 1:HCV$numberSamples){ 
  
  Param_pop_array[[set]] <- Param_pop_array[[set]][ , ,c(1:trimpt)]
  
}
### 

tic <- proc.time()

if (runSamples) {
  
  ScenparamResults <- list()

    for(m in names(CEA_mainParam)){
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
rm(Param_dfList, Param_dfListExtend, steady, CEA_main)
gc()

# parameter sets for scenarios with uptake PrEP and HIVD 
ScenparamResults_PrEP <- list()
  
ScenparamResults_HIVD <- list()
  
ScenparamResults_PrEPnHIVD <- list()

tic <- proc.time()

for(m in names(CEA_mainParam)){
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
gc()

tic <- proc.time()

for(m in names(CEA_mainParam)){
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

gc()

tic <- proc.time()

for(m in names(CEA_mainParam)){
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



