# running scenario 

# time-varying parameters 
# testing: tau_ab, tau_ag, tau_poct 
# treatment: eta 
# rm(list = ls())

library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)

project_name <- "POC_AU"

codefun_path <- paste("/Users/jjwu/Documents/Simplified-HCV-testing-model")

data_path <- paste("/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/", 
                   project_name, sep = "")
# Load useful libraries

library("readr")
library("dplyr")
library("tidyr")
library("purrr")
library("parallel")
library("pacman")
library("doMC")

Rcode <- file.path(codefun_path, "03. Code")

DataFolder <- file.path(data_path, "01. DATA/model input" )
OutputFolder <- file.path(data_path, "02. Output")

load(file.path(OutputFolder, paste0(project_name, ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_init", ".rda")))
load(file.path(OutputFolder, paste0(project_name, "cali", ".rda")))
# load(file.path(OutputFolder, paste0(project_name, "cali.rda")))


urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- FALSE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))


# % increasing for stages in cascade

treatuptake <- c(0.48, 0.48, 0.48, 0.48, 0.4)
calibration_Y <- 2015
varying_Yint <- 2015
varying_Yend <-2021
varying_tend <- 2022
varyingYpoint_int <- (varying_Yint - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_tend <- (varying_tend - calibration_Y)/POC_AU$timestep + 1

intVal <- dfList$eta[, 3, 1] 

for ( j in 1: (POC_AU$npops)){
  for ( i in 2:dim(dfList$eta)[[2]]) {

    dfList$eta[j, i, c(varyingYpoint_int:varyingYpoint_end)] <- 
      c(seq(intVal[j], treatuptake[j], length = varyingYpoint_end - varyingYpoint_int + 1)) 
    
    dfList$eta[j, i, c((varyingYpoint_end + 1): dim(dfList$eta)[[3]])] <- treatuptake[j]
  }
  
}

save(POC_AU, constants, disease_progress, fib, dfList, pop_array, 
     constantsDf, initialPops, best_estimates, best_initial_pop, 
     file = file.path(OutputFolder,
                      paste0(project_name, ".rda")))

dfList_s <- dfList
# changing the parameters 

#RNA testing 
dfList$tau_RNA[ , ,1]

intVal <- dfList$tau_RNA[, 3, 1] 
intVal_poct <- dfList$tau_poct[, 3, 1] 
endV <- c(0.5*1.36, 0.42*1.61, 0.5*1.36, 0.4*1.61, 0.4*1.61)
endV_poct <- c(0.28*2.14, 0.23*2.58, 0.1*2, 0.16*2, 0.052*2)

for ( j in 1: (POC_AU$npops)){
  for ( i in 2:dim(dfList$eta)[[2]]) {
    
    dfList_s$tau_RNA[j, i, c(varyingYpoint_end:varyingYpoint_tend)] <- 
      c(seq(intVal[j], endV[j], length = varyingYpoint_tend - varyingYpoint_end + 1)) 
    
    dfList_s$tau_RNA[j, i, c((varyingYpoint_tend + 1): dim(dfList$eta)[[3]])] <- endV[j]
    
    dfList_s$tau_poct[j, i, c(varyingYpoint_end:varyingYpoint_tend)] <- 
      c(seq(intVal_poct[j], endV_poct[j], length = varyingYpoint_tend - varyingYpoint_end + 1)) 
    
    dfList_s$tau_poct[j, i, c((varyingYpoint_tend + 1): dim(dfList$eta)[[3]])] <- endV_poct[j]
    
    
  }
  
}


tic <- proc.time()
endY <- 100
calibrateInit <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                        disease_progress,pop_array,
                        dfList, fib,
                        modelrun="UN", proj = "POC_AU", end_Y = endY)


toc <- proc.time() - tic 
toc



save(calibrateInit, 
     file = file.path(OutputFolder ,
                      paste0(project_name,"cali_timevarying" ,".rda")))


tic <- proc.time()
endY <- 100
scenario_conser <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                        disease_progress,pop_array,
                        dfList_s, fib,
                        modelrun="UN", proj = "POC_AU", end_Y = endY)


toc <- proc.time() - tic 
toc


save(scenario_conser, 
     file = file.path(OutputFolder ,
                      paste0(project_name,"scenario_conser" ,".rda")))
