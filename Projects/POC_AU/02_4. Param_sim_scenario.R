#### need to seperate the param_sim for scenarios #### 
rm(list = ls())
gc()
project_name <- "POC_AU"
options(digits = 15)

codefun_path <- paste("/Users/jjwu/Projects/Simplified-HCV-testing-model")

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
RDAFolder <- file.path(data_path, "02. Output")
OutputFolder <- file.path(codefun_path, "Projects/POC_AU/Output")
OutputFig <- file.path(codefun_path, "Projects/POC_AU/Figs")

load(file.path(RDAFolder, paste0(project_name, ".rda")))

load(file.path(RDAFolder, paste0(project_name, "param.rda")))
load(file.path(RDAFolder, paste0(project_name, "paramDflist.rda")))

source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R")) 

source(file.path(Rcode, "/Functions/plotManuscript.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R")) 

urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- TRUE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing
number_samples <- 1000

POC_AU$numberSamples <- number_samples

tic <- proc.time()
trim_pt <- 100*(1/POC_AU$timestep)
Param_estimates <- lapply(Param_estimates, function(x) x[c(1:trim_pt),]%>%as.data.frame)
param_poparray <- lapply(param_poparray , function(x) x[, , c(1:trim_pt)])
gc()
paramDflist <- lapply(paramDflist, function(x) lapply(x, function(y) y[, , c(1:trim_pt)]))
gc()

#### cost sensitivity #####
# run by each scenarios for faster result generation, rather than multiple loop
cost_types <- c("fixednvariable", "total", "DAAcost_reducquarter", "DAAcost_reduchalf")

##### scenarios ##### 
load(file.path(RDAFolder, paste0(project_name, "scenario_cascade.rda")))

sce_name <- names(scenario_cascade)
rm(scenario_cascade)

param_dfList <- list()

param_scenario <- list()
trim_pt <- 100*(1/POC_AU$timestep)

# check whether any parameter >=1 : which_ones <- which(sapply(1:1000, function(i) any(scenario_p[[i]]$tau_ab >= 1)))
# by scenarios 
load(file.path(RDAFolder, paste0(project_name,"param_scenario_",sce_name[1], ".rda"))) 

scenario_p <- lapply(scenario_p, function(x) lapply(x, function(y)y[, , c(1:trim_pt)]))
gc()
param_scenario <- list()

for(cost_type in cost_types){ 
  tic <- proc.time()
  load(file.path(RDAFolder, paste0(project_name, "param_cost_", cost_type,".rda")))
  param_scenario <- list()
  endY <- 100
  for(x in 1:1000){
    param_scenario[[x]] <- HCVMSM(POC_AU, Param_estimates[[x]], Param_Pops[[x]],
                                  Param_disease_progress[[x]], param_poparray[[x]],
                                  paramDflist[[x]], param_cascade_sc = scenario_p[[x]], 
                                  fib = Param_fib[[x]], 
                                  modelrun="UN", proj = "POC_AU", end_Y = endY, 
                                  cost = param_cost[[x]], costflow = param_cost_flow[[x]], 
                                  costflow_Neg = param_costflow_Neg[[x]], fc_sc = scenario_fc[[1]],
                                  fp = NULL)
    
  }
  toc <- proc.time() - tic
  print(paste0("Completed: ", cost_type, " | Time: ", toc))
  save(param_scenario,
       file = file.path(OutputFolder,
                        paste0(project_name, "param_sc_", sce_name[1], "_", cost_type, ".rda")))
  
  print(paste0("Saved: ", cost_type))
  
  rm(param_scenario)
  gc()
}

####Scenario 2 ####
load(file.path(RDAFolder, paste0(project_name,"param_scenario_",sce_name[2], ".rda"))) 

scenario_p <- lapply(scenario_p, function(x) lapply(x, function(y)y[, , c(1:trim_pt)]))
gc()
param_scenario <- list()

for(cost_type in cost_types){ 
  tic <- proc.time()
  load(file.path(RDAFolder, paste0(project_name, "param_cost_", cost_type,".rda")))
  param_scenario <- list()
  endY <- 100
  for(x in 1:1000){
    param_scenario[[x]] <- HCVMSM(POC_AU, Param_estimates[[x]], Param_Pops[[x]],
                                  Param_disease_progress[[x]], param_poparray[[x]],
                                  paramDflist[[x]], param_cascade_sc = scenario_p[[x]], 
                                  fib = Param_fib[[x]], 
                                  modelrun="UN", proj = "POC_AU", end_Y = endY, 
                                  cost = param_cost[[x]], costflow = param_cost_flow[[x]], 
                                  costflow_Neg = param_costflow_Neg[[x]], fc_sc = scenario_fc[[2]],
                                  fp = NULL)
    
  }
  toc <- proc.time() - tic
  print(paste0("Completed: ", cost_type, " | Time: ", toc))
  save(param_scenario,
       file = file.path(OutputFolder,
                        paste0(project_name, "param_sc_", sce_name[2], "_", cost_type, ".rda")))
  
  print(paste0("Saved: ", cost_type))
  
  rm(param_scenario)
  gc()
}

####Scenario 3 ####
load(file.path(RDAFolder, paste0(project_name,"param_scenario_",sce_name[3], ".rda"))) 

scenario_p <- lapply(scenario_p, function(x) lapply(x, function(y)y[, , c(1:trim_pt)]))
gc()
param_scenario <- list()

for(cost_type in cost_types){ 
  tic <- proc.time()
  load(file.path(RDAFolder, paste0(project_name, "param_cost_", cost_type,".rda")))
  param_scenario <- list()
  endY <- 100
  for(x in 1:1000){
    param_scenario[[x]] <- HCVMSM(POC_AU, Param_estimates[[x]], Param_Pops[[x]],
                                  Param_disease_progress[[x]], param_poparray[[x]],
                                  paramDflist[[x]], param_cascade_sc = scenario_p[[x]], 
                                  fib = Param_fib[[x]], 
                                  modelrun="UN", proj = "POC_AU", end_Y = endY, 
                                  cost = param_cost[[x]], costflow = param_cost_flow[[x]], 
                                  costflow_Neg = param_costflow_Neg[[x]], fc_sc = scenario_fc[[3]],
                                  fp = NULL)
    
  }
  toc <- proc.time() - tic
  print(paste0("Completed: ", cost_type, " | Time: ", toc))
  save(param_scenario,
       file = file.path(OutputFolder,
                        paste0(project_name, "param_sc_", sce_name[3], "_", cost_type, ".rda")))
  
  print(paste0("Saved: ", cost_type))
  
  rm(param_scenario)
  gc()
}


####Scenario 4 ####
load(file.path(RDAFolder, paste0(project_name,"param_scenario_",sce_name[4], ".rda"))) 

scenario_p <- lapply(scenario_p, function(x) lapply(x, function(y)y[, , c(1:trim_pt)]))
gc()
param_scenario <- list()

for(cost_type in cost_types){ 
  tic <- proc.time()
  load(file.path(RDAFolder, paste0(project_name, "param_cost_", cost_type,".rda")))
  param_scenario <- list()
  endY <- 100
  for(x in 1:1000){
    param_scenario[[x]] <- HCVMSM(POC_AU, Param_estimates[[x]], Param_Pops[[x]],
                                  Param_disease_progress[[x]], param_poparray[[x]],
                                  paramDflist[[x]], param_cascade_sc = scenario_p[[x]], 
                                  fib = Param_fib[[x]], 
                                  modelrun="UN", proj = "POC_AU", end_Y = endY, 
                                  cost = param_cost[[x]], costflow = param_cost_flow[[x]], 
                                  costflow_Neg = param_costflow_Neg[[x]], fc_sc = scenario_fc[[4]],
                                  fp = NULL)
    
  }
  toc <- proc.time() - tic
  print(paste0("Completed: ", cost_type, " | Time: ", toc))
  save(param_scenario,
       file = file.path(OutputFolder,
                        paste0(project_name, "param_sc_", sce_name[4], "_", cost_type, ".rda")))
  
  print(paste0("Saved: ", cost_type))
  
  rm(param_scenario)
  gc()
}
