
#===============================================================================
# S4: PrEPnHIVD
#===============================================================================

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

load(file.path(paste0(outputdt, "/simScen_PrEPHIV.rda")))

load(file.path(paste0(outputdt, "/simScen_param_PrEPnHIVD.rda")))

rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile) 

tic <- proc.time()

Outcome_Scen_PrEPnHIVD <- list()

ScenResults_PrEPHIV$PrEPnHIVD <- ScenResults_PrEPHIV$PrEPnHIVD[names(ScenResults_PrEPHIV$PrEPnHIVD)!="DBS"]

ScenparamResults_PrEPnHIVD <- ScenparamResults_PrEPnHIVD[names(ScenparamResults_PrEPnHIVD) != "DBS"] 

Outcome_Scen_PrEPnHIVD[["epi"]] <- Res_summary_epi(HCV, 
                                              ScenResults_PrEPHIV$PrEPnHIVD, 
                                              ScenparamResults_PrEPnHIVD, 
                                              endY = 100, base_case = NULL)


Outcome_Scen_PrEPnHIVD[["costqaly"]] <- lapply(names(ScenResults_PrEPHIV$PrEPnHIVD), 
                                      function(i) Res_summary_costqaly(HCV, 
                                                                       ScenResults_PrEPHIV$PrEPnHIVD[[i]], 
                                                                       ScenparamResults_PrEPnHIVD[[i]], 
                                                                       endY = 100)
)

names(Outcome_Scen_PrEPnHIVD[["costqaly"]]) <- names(ScenResults_PrEPHIV$PrEPnHIVD)

toc <- proc.time() - tic

toc

save(Outcome_Scen_PrEPnHIVD,
     file = file.path(outputdt, "Outcome_Scen_PrEPnHIVD.rda"))

rm(Outcome_Scen_PrEPnHIVD, ScenparamResults_PrEPnHIVD)

gc()

