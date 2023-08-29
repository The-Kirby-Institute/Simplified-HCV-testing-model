#===============================================================================
# S3: HIVD
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

load(file.path(paste0(outputdt, "/simScen_param_HIVD.rda")))

rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile) 

tic <- proc.time()

Outcome_Scen_HIVD <- list()

ScenResults_PrEPHIV$HIVD<- ScenResults_PrEPHIV$HIVD[names(ScenResults_PrEPHIV$HIVD)!="DBS"]

ScenparamResults_HIVD <- ScenparamResults_HIVD[names(ScenparamResults_HIVD) != "DBS"] 

Outcome_Scen_HIVD[["epi"]] <- Res_summary_epi(HCV, 
                                         ScenResults_PrEPHIV$HIVD, 
                                         ScenparamResults_HIVD, 
                                         endY = 100, base_case = NULL)


Outcome_Scen_HIVD[["costqaly"]] <- lapply(names(ScenResults_PrEPHIV$HIVD), 
                                 function(i) Res_summary_costqaly(HCV, 
                                                                  ScenResults_PrEPHIV$HIVD[[i]], 
                                                                  ScenparamResults_HIVD[[i]], 
                                                                  endY = 100)
)

names(Outcome_Scen_HIVD[["costqaly"]]) <- names(ScenResults_PrEPHIV$HIVD)

toc <- proc.time() - tic

toc

save(Outcome_Scen_HIVD,
     file = file.path(outputdt, "Outcome_Scen_HIVD.rda"))

rm(Outcome_Scen_HIVD, ScenparamResults_HIVD)

gc()
