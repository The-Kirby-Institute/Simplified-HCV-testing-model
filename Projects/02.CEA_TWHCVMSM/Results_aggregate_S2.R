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

codep <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/TWHCV-model"

com_codeP <- "/Users/jjwu/Documents/Simplified-HCV-testing-model/03. Code"

epidatapath <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Taiwan-MSM-HCV-model"

project_codep <- "/Users/jjwu/Documents/Simplified-HCV-testing-model/Projects/02.CEA_TWHCVMSM"

#file path of "TWHCV-model" project
Rcode <- file.path(codep, '03. Code')

scenariopath <- file.path(epidatapath, '01. DATA/model input/Scenarios')

# save and cost data folder 
dtp <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/02.CEA_TWHCVMSM"

DataFolder <- file.path(dtp, "01. DATA/model input")

outputdt <- file.path(dtp, "02. Output/RDA")

outputfig <- file.path(dtp, "02. Output/Figs")

# source 

source(file.path(com_codeP, "Functions/plotFunctions.R"))

source(file.path(com_codeP, "Functions/plotOptions.R"))

source(file.path(com_codeP, "Functions/Scenarios_wrangling.R"))

source(file.path(com_codeP, "Functions/plotManuscript.R"))

source(file.path(project_codep, "AggregateRes.R"))


load(file.path(paste0(outputdt, "/simScen_PrEPHIV.rda")))

load(file.path(paste0(outputdt, "/simScen_param_PrEP.rda")))

projectFile <- file.path(epidatapath , paste0("HCVModel",".rda"))

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

#### modify #### 

x <- lapply(names(ScenResults_PrEPHIV$PrEP), 
            function(i) Res_summary_costqaly(HCV, 
                                             ScenResults_PrEPHIV$PrEP[[i]], 
                                             ScenparamResults_PrEP[[i]], 
                                             endY = 100))

names(x) <- names(ScenResults_PrEPHIV$PrEP)

rm(ScenResults_PrEPHIV$PrEP, ScenparamResults_PrEP)

gc()

load(file.path(paste0(outputdt, "/Outcome_Scen_PrEP.rda")))

Outcome_Scen_PrEP$costqaly <- x

save(Outcome_Scen_PrEP,
     file = file.path(outputdt, "Outcome_Scen_PrEP.rda"))

rm(Outcome_Scen_PrEP, ScenparamResults_PrEP)

gc()
