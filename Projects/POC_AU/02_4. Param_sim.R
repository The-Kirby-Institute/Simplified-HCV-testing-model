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

load(file.path(OutputFolder, paste0(project_name, ".rda")))

load(file.path(OutputFolder, paste0(project_name, "param.rda")))
load(file.path(OutputFolder, paste0(project_name, "paramDflist.rda")))
load(file.path(OutputFolder, paste0(project_name, "param_cost.rda")))
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





param_sq <- list()
param_dfList <- lapply(paramDflist[[1]], function(x) x*0)

names(param_dfList) <- names(paramDflist[[1]])

fc <- matrix(0, ncol = dim(paramDflist[[1]]$eta)[3], nrow = POC_AU$npops)
endY <- 100
tic <- proc.time()
trim_pt <- 100*(1/POC_AU$timestep)
Param_estimates <- lapply(Param_estimates, function(x) x[c(1:trim_pt),]%>%as.data.frame)
param_poparray <- lapply(param_poparray , function(x) x[, , c(1:trim_pt)])

paramDflist <- lapply(paramDflist, function(x) lapply(x, function(y) y[, , c(1:trim_pt)]))
gc()

for(x in 1:1000){
  param_sq[[x]] <- HCVMSM(POC_AU, Param_estimates[[x]], Param_Pops[[x]],
                          Param_disease_progress[[x]], param_poparray[[x]],
                          paramDflist[[x]], param_cascade_sc = param_dfList, 
                          fib = Param_fib[[x]], 
                          modelrun="UN", proj = "POC_AU", end_Y = endY, 
                          cost = param_cost[[x]], costflow = param_cost_flow[[x]], 
                          costflow_Neg = param_costflow_Neg[[x]], fc_sc = fc,
                          fp = NULL)
  

  
  
}

toc <- proc.time() - tic

save(param_sq,
     file = file.path(OutputFolder,
                      paste0(project_name, "param_simulation", ".rda")))

rm(param_sq) 
gc()
##### scenarios ##### 
load(file.path(OutputFolder, paste0(project_name, "scenario_cascade.rda")))

sce_name <- names(scenario_cascade)[!names(scenario_cascade)%in%
                                                             c("dfList_NPexp_B", "dfList_NPexp_C")] 

rm(scenario_cascade)

param_dfList <- list()

param_scenario <- list()
trim_pt <- 100*(1/POC_AU$timestep)
for(n in sce_name){ 
  load(file.path(OutputFolder, paste0(project_name,"param_scenario_",n, ".rda"))) 
  
  trim_pt <- 100*(1/POC_AU$timestep)
  scenario_p <- lapply(scenario_p, function(x) lapply(x, function(y)y[, , c(1:trim_pt)]))
  gc()
  tic <- proc.time()
  param_scenario <- list()
  for(x in 1:1000){
    param_scenario[[x]] <- HCVMSM(POC_AU, Param_estimates[[x]], Param_Pops[[x]],
                                  Param_disease_progress[[x]], param_poparray[[x]],
                                  paramDflist[[x]], param_cascade_sc = scenario_p[[x]], 
                                  fib = Param_fib[[x]], 
                                  modelrun="UN", proj = "POC_AU", end_Y = endY, 
                                  cost = param_cost[[x]], costflow = param_cost_flow[[x]], 
                                  costflow_Neg = param_costflow_Neg[[x]], fc_sc = scenario_fc[[n]],
                                  fp = NULL)
    
  }
  
  toc <- proc.time() - tic
  
  save(param_scenario,
       file = file.path(OutputFolder,
                        paste0(project_name, "param_sc_", n, ".rda")))
  
  rm(scenario_p, param_scenario)
  gc()
}





