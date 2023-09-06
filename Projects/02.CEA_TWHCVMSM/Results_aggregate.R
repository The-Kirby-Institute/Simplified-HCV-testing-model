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

files <- list.files(path = paste0(outputdt, "/", sep = ""),
                    pattern = '^simBase.*\\.rda')
# glue the files name and path 
files <- lapply(files, function(f) paste0(file.path(paste0(outputdt, "/" ,f, sep =""))))

# load to gloval environment 
sapply(files, function(f) {load(f, envir=globalenv())})

# Rda file path 
# load the .rda file of base estimate 
projectFile <- file.path(epidatapath , paste0("HCVModel",".rda"))

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

gc()

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
str(SQ)
# test

################################################################################
 tic <- proc.time()

Outcome_Base_cost <- list() 

Outcome_Base_cost <- lapply(names(SQ), function(i) 
  Res_summary_costqaly(HCV, 
                       SQ[[i]], 
                       SQ_param[[i]], 
                       endY = 100)
  )
names(Outcome_Base_cost) <- names(SQ)


# reorganizing list 

Outcome_Base_PrEPsame <- list() 



save(Outcome_Base_epi,
     file = file.path(outputdt, "Outcome_Base_epi.rda"))

save(Outcome_Base_cost,
     file = file.path(outputdt, "Outcome_Base_cost.rda"))

rm(Outcome_Base_epi, Outcome_Base_cost, SQ, SQ_param)



load(file.path(paste0(outputdt, "/Outcome_Base_PrEPsame.rda")))

Outcome_Base_PrEPsame$costqaly <- Outcome_Base_cost$Base
save(Outcome_Base_PrEPsame,
     file = file.path(outputdt, "Outcome_Base_PrEPsame.rda"))
rm(Outcome_Base_PrEPsame)


load(file.path(paste0(outputdt, "/Outcome_Base_PrEP.rda")))

Outcome_Base_PrEP$costqaly <- Outcome_Base_cost$PrEP
  save(Outcome_Base_PrEP,
       file = file.path(outputdt, "Outcome_Base_PrEP.rda"))
rm(Outcome_Base_PrEP)

load(file.path(paste0(outputdt, "/Outcome_Base_HIVD.rda")))

Outcome_Base_HIVD$costqaly <- Outcome_Base_cost$HIVD
save(Outcome_Base_HIVD,
     file = file.path(outputdt, "Outcome_Base_HIVD.rda"))
rm(Outcome_Base_HIVD)

load(file.path(paste0(outputdt, "/Outcome_Base_PrEPnHIVD.rda")))

Outcome_Base_PrEPnHIVD$costqaly <- Outcome_Base_cost$PrEPnHIVD
save(Outcome_Base_PrEPnHIVD,
     file = file.path(outputdt, "Outcome_Base_PrEPnHIVD.rda"))
rm(Outcome_Base_PrEPnHIVD)
#### load outcome file and append the new costqaly list 



