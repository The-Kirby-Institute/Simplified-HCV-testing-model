##### 
# Things to do 
# saperate scenarios to saperate script.
# restart before running the next script to aviod R crush
# consider do these from the terminal 
#==============================================================================
# status quo in each PrEP coverage scenario 
rm(list = ls()) 

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

files <- list.files(path = paste0(outputdt, "/", sep = ""),
                    pattern = '^simBase.*\\.rda')
# glue the files name and path 
files <- lapply(files, function(f) paste0(file.path(paste0(outputdt, "/" ,f, sep =""))))

# load to gloval environment 
sapply(files, function(f) {load(f)})

# Rda file path 
# load the .rda file of base estimate 
rdapath <- file.path(here()%>%dirname(), "Taiwan-MSM-HCV-model")

projectFile <- file.path(rdapath , paste0("HCVModel",".rda"))

projectVars <- load(projectFile) 

# combine simulations in a list 

SQ <- list(Base = bestResults, 
           PrEP = bestResults_PrEPHIV$PrEP,
           HIVD = bestResults_PrEPHIV$HIVD,
           PrEPnHIVD = bestResults_PrEPHIV$PrEPnHIVD)

rm(bestResults, bestResults_PrEPHIV)


SQ_param <- list(Base = paramResults, 
                 PrEP = paramResults_PrEP,
                 HIVD = paramResults_HIVD,
                 PrEPnHIVD = paramResults_PrEPnHIVD)

rm(paramResults, paramResults_HIVD, paramResults_PrEP, paramResults_PrEPnHIVD)


# generate the results dataframes 
tic <- proc.time()

Outcome_Base_epi <- list()

Outcome_Base_epi <- lapply(names(SQ), function(i) 
  Res_summary_epi(HCV, 
                  SQ[[i]], 
                  SQ_param[[i]], 
                  endY = 100, base_case = "y"))

names(Outcome_Base_epi) <- names(SQ)

toc <- proc.time() - tic 

toc 


 tic <- proc.time()

Outcome_Base_cost <- list() 

Outcome_Base_cost <- lapply(names(SQ), function(i) 
  Res_summary_costqaly(HCV, 
                       SQ[[i]], 
                       SQ_param[[i]], 
                       endY = 100)
  )
names(Outcome_Base_cost) <- names(SQ)

save(Outcome_Base_epi,
     file = file.path(outputdt, "Outcome_Base_epi.rda"))

save(Outcome_Base_cost,
     file = file.path(outputdt, "Outcome_Base_cost.rda"))

rm(Outcome_Base_epi, Outcome_Base_cost, SQ, SQ_param)



######### 
# S1: HIV services coverage remain unchanged 

load(file.path(paste0(outputdt, "/simScen.rda")))

tic <- proc.time()

Outcome_Scen_epi_Base <- list()

  
Outcome_Scen_epi_Base <- Res_summary_epi(HCV, 
                                       ScenResults, 
                                       ScenparamResults, 
                                       endY = 100, base_case = NULL)
  
Outcome_Scen_cost_Base <- list()

Outcome_Scen_cost_Base <- lapply(names(ScenResults), 
                                 function(i) Res_summary_costqaly(HCV, 
                                                                  ScenResults[[i]], 
                                                                  ScenparamResults[[i]], 
                                                                  endY = 100)
                                 )

names(Outcome_Scen_cost_Base) <- names(ScenResults)

toc <- proc.time() - tic

toc

save(Outcome_Scen_epi_Base,
     file = file.path(outputdt, "Outcome_Scen_epi_Base.rda"))

save(Outcome_Scen_cost_Base,
     file = file.path(outputdt, "Outcome_Scen_cost_Base.rda"))

rm(Outcome_Scen_epi_Base, Outcome_Scen_cost_Base)

#===============================================================================
# S2-S4 
# the rda files are highly memory demand 
# requires to restart for each scenarios 

rm(list = ls()) 

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

Outcome_Scen_epi_PrEP <- list()

Outcome_Scen_epi_PrEP <- Res_summary_epi(HCV, 
                                         ScenResults_PrEPHIV$PrEP, 
                                         ScenparamResults_PrEP, 
                                         endY = 100, base_case = NULL)

Outcome_Scen_cost_PrEP <- list()

Outcome_Scen_cost_PrEP <- lapply(names(ScenResults_PrEPHIV$PrEP), 
                                 function(i) Res_summary_costqaly(HCV, 
                                                                  ScenResults_PrEPHIV$PrEP[[i]], 
                                                                  ScenparamResults_PrEP[[i]], 
                                                                  endY = 100)
)

names(Outcome_Scen_cost_PrEP) <- names(ScenResults_PrEPHIV$PrEP)

toc <- proc.time() - tic

toc

save(Outcome_Scen_epi_PrEP,
     file = file.path(outputdt, "Outcome_Scen_epi_PrEP.rda"))

save(Outcome_Scen_cost_PrEP,
     file = file.path(outputdt, "Outcome_Scen_cost_PrEP.rda"))

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

Outcome_Scen_epi_HIVD <- list()

Outcome_Scen_epi_HIVD <- Res_summary_epi(HCV, 
                                         ScenResults_PrEPHIV$HIVD, 
                                         ScenparamResults_HIVD, 
                                         endY = 100, base_case = NULL)

Outcome_Scen_cost_HIVD <- list()

Outcome_Scen_cost_HIVD <- lapply(names(ScenResults_PrEPHIV$HIVD), 
                                 function(i) Res_summary_costqaly(HCV, 
                                                                  ScenResults_PrEPHIV$HIVD[[i]], 
                                                                  ScenparamResults_HIVD[[i]], 
                                                                  endY = 100)
)

names(Outcome_Scen_cost_HIVD) <- names(ScenResults_PrEPHIV$HIVD)

toc <- proc.time() - tic

toc

save(Outcome_Scen_epi_HIVD,
     file = file.path(outputdt, "Outcome_Scen_epi_HIVD.rda"))

save(Outcome_Scen_cost_HIVD,
     file = file.path(outputdt, "Outcome_Scen_cost_HIVD.rda"))


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

Outcome_Scen_epi_PrEPnHIVD <- list()

Outcome_Scen_epi_PrEPnHIVD <- Res_summary_epi(HCV, 
                                         ScenResults_PrEPHIV$PrEPnHIVD, 
                                         ScenparamResults_PrEPnHIVD, 
                                         endY = 100, base_case = NULL)

Outcome_Scen_cost_PrEPnHIVD <- list()

Outcome_Scen_cost_PrEPnHIVD <- lapply(names(ScenResults_PrEPHIV$PrEPnHIVD), 
                                 function(i) Res_summary_costqaly(HCV, 
                                                                  ScenResults_PrEPHIV$PrEPnHIVD[[i]], 
                                                                  ScenparamResults_PrEPnHIVD[[i]], 
                                                                  endY = 100)
)

names(Outcome_Scen_cost_PrEPnHIVD) <- names(ScenResults_PrEPHIV$PrEPnHIVD)

toc <- proc.time() - tic

toc

save(Outcome_Scen_epi_PrEPnHIVD,
     file = file.path(outputdt, "Outcome_Scen_epi_PrEPnHIVD.rda"))

save(Outcome_Scen_cost_PrEPnHIVD,
     file = file.path(outputdt, "Outcome_Scen_cost_PrEPnHIVD.rda"))




