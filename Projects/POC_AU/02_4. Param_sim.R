rm(list = ls())

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
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
load(file.path(OutputFolder, paste0(project_name, "param.rda")))
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
start_lower <- 0.75
start_upper <- 1.25
end_lower <- 0.75
end_upper <- 1.25


POC_AU$numberSamples <- number_samples


# cost dt
files <- list.files(path = paste0(DataFolder, 
                                  "/cost/", sep =  ""), pattern = '*.csv')


costdfList <- lapply(files, function(f) {
  
  df <- read.csv(file.path(paste0(DataFolder, "/cost/", f, sep = "")), header = TRUE)
  
  df <- df[, -1]
  
  df <- df%>%as_tibble()
  
  df <- as.matrix(df, nrow = npops, ncol = length(.) + 1)
  
})

names(costdfList) <- c(gsub("^|.csv", "", files)) # ^: from beginning, \ end before .csv


cost_state <- costdfList$state
costflow <- list()
costflow[[1]] <- costdfList$costFlow
costflow[[2]] <- costdfList$costFlow_POCRNA

costflow_Neg <- list()
costflow_Neg[[1]] <- costdfList$costFlow_NEG
costflow_Neg[[2]] <- costdfList$`costFlow_POCRNA _NEG`




param_sq <- list()
param_dfList <- lapply(dfList, function(x) x*0)
names(param_dfList) <- names(dfList)

fc <- matrix(0, ncol = dim(dfList$eta)[3], nrow = POC_AU$npops)
endY <- 100
tic <- proc.time()
trim_pt <- 100*(1/POC_AU$timestep)
Param_estimates <- lapply(Param_estimates, function(x) x[c(1:trim_pt),])
param_poparray <- lapply(param_poparray , function(x) x[, , c(1:trim_pt)])
paramDflist <- lapply(paramDflist, function(x) lapply(x, function(y) y[, , c(1:trim_pt)]))
gc()
for(i in 1:1000){
  param_sq[[i]] <- HCVMSM(POC_AU, Param_estimates[[i]], best_est_pop,
                          Param_disease_progress[[i]], param_poparray[[i]],
                          paramDflist[[i]], param_cascade_sc = param_dfList , 
                          fib = Param_fib[[i]], 
                          modelrun="UN", proj = "POC_AU", end_Y = endY, 
                          cost = costdfList, costflow = costflow, 
                          costflow_Neg = costflow_Neg, fc_sc = fc,
                          fp = NULL)
  

  
  
}

toc <- proc.time() - tic

save(param_sq,
     file = file.path(OutputFolder,
                      paste0(project_name, "param_simulation", ".rda")))


