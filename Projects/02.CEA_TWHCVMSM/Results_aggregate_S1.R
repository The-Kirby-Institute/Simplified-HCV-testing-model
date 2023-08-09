#### Result_aggregate
####
# S1: PrEP coverage remained unchanged
####
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

load(file.path(paste0(outputdt, "/simScen.rda")))

rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile) 

tic <- proc.time()

Outcome_Scen_PrEPsame <- list()

ScenResults <- ScenResults[names(ScenResults)!="DBS"]

ScenparamResults <- ScenparamResults[names(ScenparamResults) != "DBS"] 

Outcome_Scen_PrEPsame[["epi"]] <- Res_summary_epi(HCV, 
                                         ScenResults, 
                                         ScenparamResults, 
                                         endY = 100, base_case = NULL)

Outcome_Scen_PrEPsame[["costqaly"]] <- lapply(names(ScenResults), 
                                 function(i) Res_summary_costqaly(HCV, 
                                                                  ScenResults[[i]], 
                                                                  ScenparamResults[[i]], 
                                                                  endY = 100))

names(Outcome_Scen_PrEPsame[["costqaly"]]) <- names(ScenResults)

toc <- proc.time() - tic

toc

save(Outcome_Scen_PrEPsame,
     file = file.path(outputdt, "Outcome_Scen_PrEPsame.rda"))


rm(Outcome_Scen_Base, ScenResults, ScenparamResults)
gc()
