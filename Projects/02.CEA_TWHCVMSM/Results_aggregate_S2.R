#===============================================================================
# S2 
# the rda files are highly memory demand 
# requires to restart for each scenarios 
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

load(file.path(paste0(outputdt, "/simScen_param_PrEP.rda")))

rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile) 

# PrEP
tic <- proc.time()

Outcome_Scen_PrEP <- list()

ScenResults_PrEPHIV$PrEP <- ScenResults_PrEPHIV$PrEP[names(ScenResults_PrEPHIV$PrEP)!="DBS"]

ScenparamResults_PrEP <- ScenparamResults_PrEP[names(ScenparamResults_PrEP) != "DBS"] 

Outcome_Scen_PrEP[["epi"]] <- Res_summary_epi(HCV, 
                                         ScenResults_PrEPHIV$PrEP, 
                                         ScenparamResults_PrEP, 
                                         endY = 100, base_case = NULL)

Outcome_Scen_PrEP[["costqaly"]] <- lapply(names(ScenResults_PrEPHIV$PrEP), 
                                 function(i) Res_summary_costqaly(HCV, 
                                                                  ScenResults_PrEPHIV$PrEP[[i]], 
                                                                  ScenparamResults_PrEP[[i]], 
                                                                  endY = 100)
)

names(Outcome_Scen_PrEP[["costqaly"]]) <- names(ScenResults_PrEPHIV$PrEP)

toc <- proc.time() - tic

toc

save(Outcome_Scen_PrEP,
     file = file.path(outputdt, "Outcome_Scen_PrEP.rda"))

rm(Outcome_Scen_PrEP, ScenparamResults_PrEP)

gc()
