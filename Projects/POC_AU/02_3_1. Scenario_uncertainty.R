# uncertainty around HCV cascade in scenarios 
rm(list = ls())
gc()
project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Documents/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
# Load useful libraries
library("lhs")
library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")
library("ggplot2")
library("viridis") 
library("openxlsx")
library(gt)
library(dplyr)
Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")
OutputFig <- file.path(OutputFolder, "Figs")


source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R")) 

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 

load(file.path(OutputFolder, paste0(project_name, "lhs_sampling.rda")))
load(file.path(OutputFolder, paste0(project_name, ".rda")))
gc() 


# function to generate parameters set for HCV cascade in scenarios 
# lota and cured did not involve in sampling, the effect will sampling around SVR
paramset_scenario <- function(paramDflist, lhs, dflist_LL, dflist_UU){ 
  for(i in 1:POC_AU$numberSamples){
    paramDflist[[i]][["cured"]] <- (dflist_LL[["cured"]] + dflist_UU[["cured"]])/2
  
  }
  for(i in 1:POC_AU$numberSamples){
    paramDflist[[i]][["eta"]] <- lhs[i,"poparray"]*(dflist_UU[["eta"]] - dflist_LL[["eta"]]) + dflist_LL[["eta"]] 
    paramDflist[[i]][["eta"]][paramDflist[[i]][["eta"]]>1] <- 1                               
  } 
  
  for(i in 1:POC_AU$numberSamples){
    paramDflist[[i]][["lota"]] <- (dflist_LL[["lota"]] + dflist_UU[["lota"]])/2
  }
  for(i in 1:POC_AU$numberSamples){
    paramDflist[[i]][["rho"]] <- lhs[i,"poparray"]*(dflist_UU[["rho"]] - dflist_LL[["rho"]]) + dflist_LL[["rho"]] 
    paramDflist[[i]][["rho"]][paramDflist[[i]][["rho"]]>1] <- 1  
  }
  for(i in 1:POC_AU$numberSamples){
    paramDflist[[i]][["tau_ab"]] <- lhs[i,"poparray"]*(dflist_UU[["tau_ab"]] - dflist_LL[["tau_ab"]]) + dflist_LL[["tau_ab"]] 
    paramDflist[[i]][["tau_ab"]][paramDflist[[i]][["tau_ab"]]>1] <- 1  
  }
  for(i in 1:POC_AU$numberSamples){  
    paramDflist[[i]][["tau_poct"]] <- lhs[i,"poparray"]*(dflist_UU[["tau_poct"]] - dflist_LL[["tau_poct"]]) + dflist_LL[["tau_poct"]] 
    paramDflist[[i]][["tau_poct"]][paramDflist[[i]][["tau_poct"]]>1] <- 1  
  }
  for(i in 1:POC_AU$numberSamples){   
    paramDflist[[i]][["tau_RNA"]] <- lhs[i,"poparray"]*(dflist_UU[["tau_RNA"]] - dflist_LL[["tau_RNA"]]) + dflist_LL[["tau_RNA"]] 
    paramDflist[[i]][["tau_RNA"]][paramDflist[[i]][["tau_RNA"]]>1] <- 1
  }
  return(paramDflist)
} 

attach(file.path(OutputFolder, paste0(project_name, "scenario_cascade.rda")))

sce_name <- names(scenario_cascade)[!names(scenario_cascade)%in%
                                      c("dfList_NPexp_B", "dfList_NPexp_C")] 

start_lower <- 0.75

start_upper <- 1.25

param_scenario <- list()
scenario_param <- list()
scenario_p <- list()
for(i in sce_name){ 
  param_scenario[[i]] <- scenario_cascade[[i]]
  
  scenario_dfList_LL <- lapply(param_scenario[[i]], function(x) x*start_lower)
  scenario_dfList_UU <- lapply(param_scenario[[i]], function(x) x*start_upper)
  
  scenario_param <- rep(list(param_scenario[[i]]), POC_AU$numberSamples) 
  
  scenario_p <- paramset_scenario(
    paramDflist = scenario_param,
    lhs = lhs_samples, 
    dflist_LL = scenario_dfList_LL, 
    dflist_UU = scenario_dfList_UU)
  
  
  save(scenario_p, file = file.path(OutputFolder,
                                                paste0(project_name, "param_scenario_",i, ".rda")))
  
  rm(scenario_p)
  gc()

  }


