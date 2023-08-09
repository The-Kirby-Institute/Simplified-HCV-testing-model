##### 
# Things to do 
# saperate scenarios to saperate script.
# restart before running the next script to aviod R crush
# consider do these from the terminal 
#==============================================================================
# status quo in each PrEP coverage scenario 
rm(list = ls()) 
gc()
library(here)
here()

# Load useful libraries
library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("gridExtra")
library(doParallel)
library("ggthemes")
registerDoParallel(cores = detectCores() - 1)
# Setup directories after setting working directory to source file 
# directory 

#file path of "TWHCV-model" project
codepath <- file.path(here() %>% dirname(), 'TWHCV-model/03. Code/Functions')

DataFolder <- file.path(here(), "01. DATA/model input")

ResultsFolder <- file.path(here(), "02. Output/Results")

outputdt <- here("02. Output/RDA")

outputfig <- here("02. Output/Figs")

# source 
source(file.path(codepath, "HCV_model.R"))

source(file.path(codepath, "plotFunctions.R"))

source(file.path(codepath, "plotOptions.R"))

source(file.path(codepath, "plotManuscript.R")) 

source(file.path(here(),"AggregateRes.R"))

load(file.path(here("02. Output/RDA"), "simBase.rda"))



# Rda file path 
# load the .rda file of base estimate 
rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile) 


# combine simulations in a list 

# generate the results dataframes 
tic <- proc.time()

Outcome_Base_PrEPsame <- list()

Outcome_Base_PrEPsame[["epi"]] <- Res_summary_epi(HCV, 
                                             bestResults, 
                                             paramResults, 
                                             endY = 100, base_case = "y")

Outcome_Base_PrEPsame[["costqaly"]] <- Res_summary_costqaly(HCV,
                                                   bestResults, 
                                                   paramResults, 
                                                   endY = 100)

toc <- proc.time() - tic 
toc

save(Outcome_Base_PrEPsame,
     file = file.path(outputdt, "Outcome_Base_PrEPsame.rda"))

rm(Outcome_Base_PrEPsame, bestResults, paramResults)

gc()
#### SQ: PrEP #### 
load(file.path(here("02. Output/RDA"), "simBase_PrEPHIV.rda"))

load(file.path(here("02. Output/RDA"), "simBase_param_PrEP.rda"))

tic <- proc.time()

Outcome_Base_PrEP <- list()

Outcome_Base_PrEP[["epi"]] <- Res_summary_epi(HCV, 
                                              bestResults_PrEPHIV[["PrEP"]], 
                                              paramResults_PrEP, 
                                                  endY = 100, base_case = "y")

Outcome_Base_PrEP[["costqaly"]] <- Res_summary_costqaly(HCV,
                                                        bestResults_PrEPHIV[["PrEP"]], 
                                                        paramResults_PrEP, 
                                                        endY = 100)


toc <- proc.time() - tic 

toc

save(Outcome_Base_PrEP,
     file = file.path(outputdt, "Outcome_Base_PrEP.rda"))

rm(Outcome_Base_PrEP, paramResults_PrEP)

gc()
#### SQ: HIVD ####
load(file.path(here("02. Output/RDA"), "simBase_param_HIVD.rda"))

tic <- proc.time()

Outcome_Base_HIVD <- list()

Outcome_Base_HIVD[["epi"]] <- Res_summary_epi(HCV, 
                                              bestResults_PrEPHIV[["HIVD"]], 
                                              paramResults_HIVD, 
                                              endY = 100, base_case = "y")

Outcome_Base_HIVD[["costqaly"]] <- Res_summary_costqaly(HCV,
                                                        bestResults_PrEPHIV[["HIVD"]], 
                                                        paramResults_HIVD, 
                                                        endY = 100)


toc <- proc.time() - tic 

toc

save(Outcome_Base_HIVD,
     file = file.path(outputdt, "Outcome_Base_HIVD.rda"))


rm(Outcome_Base_HIVD, paramResults_HIVD)

gc()
#### SQ: PrEPnHIVD ####
load(file.path(here("02. Output/RDA"), "simBase_param_PrEPnHIVD.rda"))

tic <- proc.time()

Outcome_Base_PrEPnHIVD <- list()

Outcome_Base_PrEPnHIVD[["epi"]] <- Res_summary_epi(HCV, 
                                              bestResults_PrEPHIV[["PrEPnHIVD"]], 
                                              paramResults_PrEPnHIVD, 
                                              endY = 100, base_case = "y")

Outcome_Base_PrEPnHIVD[["costqaly"]] <- Res_summary_costqaly(HCV,
                                                        bestResults_PrEPHIV[["PrEPnHIVD"]], 
                                                        paramResults_PrEPnHIVD, 
                                                        endY = 100)


toc <- proc.time() - tic 

toc

save(Outcome_Base_PrEPnHIVD,
     file = file.path(outputdt, "Outcome_Base_PrEPnHIVD.rda"))


rm(Outcome_Base_PrEPnHIVD, paramResults_PrEPnHIVD)

gc()
