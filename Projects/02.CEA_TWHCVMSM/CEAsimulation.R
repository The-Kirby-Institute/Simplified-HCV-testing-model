# runnuing scenarios best fit 
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

tic <- proc.time()

# Run model on best estimates

bestResults <- HCVMSM(HCV,best_estimatesOff, best_est_pop,
                      disease_progress,pop_array, dfList, fib,
                      end_Y = simY, modelrun="UN", cost = costTrim, 
                      costflow = costTrim$flow)

rm(steady)


# for scenarios related to uptake PrEP and HIV diagnosis rate
bestResults_PrEPHIV <- list()

for(i in names(pop_arrayPrEPHIV)){
  
  bestResults_PrEPHIV[[i]] <- HCVMSM(HCV,best_estimatesOff, best_est_pop,
                        disease_progress,pop_arrayPrEPHIV[[i]], dfList, fib,
                          end_Y = simY, modelrun="UN", cost = costTrim, 
                        costflow = costTrim$flow)
  
  
  
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
  
  paramResults <- list()
  
  for (set in 1:HCV$numberSamples){
  
        paramResults[[set]] <- HCVMSM(HCV,
                                  Param_estimatesOff[[set]], 
                                  Param_Pops[[set]], 
                                  Param_disease_progress[[set]], 
                                  Param_pop_array[[set]], 
                                  Param_dfListExtend[[set]],
                                  Param_fib[[set]],
                                  end_Y = simY,modelrun="UN",
                                  cost = Param_costList[[set]], 
                                  costflow = Param_costList[[set]]$flow)
        }
  }

toc <- proc.time() - tic

toc

save(bestResults, paramResults,
     file = file.path(outputdt, "simBase.rda"))

save(bestResults_PrEPHIV,
     file = file.path(outputdt, "simBase_PrEPHIV.rda"))
  
rm(bestResults, paramResults, bestResults_PrEPHIV)


tic <- proc.time()

if (runSamples) {
  
  paramResults_PrEP <- list()
  
  paramResults_HIVD <- list()
  
  paramResults_PrEPnHIVD <- list()
  
  for (set in 1:HCV$numberSamples){
    
    paramResults_PrEP[[set]] <- HCVMSM(HCV,
                                  Param_estimatesOff[[set]], 
                                  Param_Pops[[set]], 
                                  Param_disease_progress[[set]], 
                                  Param_pop_arrayPrEPHIV[["PrEP"]][[set]], 
                                  Param_dfListExtend[[set]],
                                  Param_fib[[set]],
                                  end_Y = simY,modelrun="UN",
                                  cost = Param_costList[[set]], 
                                  costflow = Param_costList[[set]]$flow)
    
    paramResults_HIVD[[set]] <- HCVMSM(HCV,
                                       Param_estimatesOff[[set]], 
                                       Param_Pops[[set]], 
                                       Param_disease_progress[[set]], 
                                       Param_pop_arrayPrEPHIV[["HIVD"]][[set]], 
                                       Param_dfListExtend[[set]],
                                       Param_fib[[set]],
                                       end_Y = simY, modelrun="UN",
                                       cost = Param_costList[[set]], 
                                       costflow = Param_costList[[set]]$flow)
    
    paramResults_PrEPnHIVD[[set]] <- HCVMSM(HCV,
                                       Param_estimatesOff[[set]], 
                                       Param_Pops[[set]], 
                                       Param_disease_progress[[set]], 
                                       Param_pop_arrayPrEPHIV[["PrEPnHIVD"]][[set]], 
                                       Param_dfListExtend[[set]],
                                       Param_fib[[set]],
                                       end_Y = simY,modelrun="UN",
                                       cost = Param_costList[[set]], 
                                       costflow = Param_costList[[set]]$flow)
    
  }
}

toc <- proc.time() - tic

toc

save(paramResults_PrEP,
     file = file.path(outputdt, "simBase_param_PrEP.rda"))

save(paramResults_HIVD,
     file = file.path(outputdt, "simBase_param_HIVD.rda"))

save(paramResults_PrEPnHIVD,
     file = file.path(outputdt, "simBase_param_PrEPnHIVD.rda"))



