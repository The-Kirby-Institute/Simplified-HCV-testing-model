# this script is for calculating the effect of national program 
# 1. We match the gap of treatment: SQ/national program + SQ
# 2. We apply a fixed of national progam coverage to match the gap 

rm(list = ls())
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
load(file.path(OutputFolder, paste0(project_name, "cali.rda")))
load(file.path(OutputFolder, paste0(project_name, "cali_timev.rda")))
source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotManuscript.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))


# 3. testing and treat parameters for national program 
RR_tauRNA_C <- 0.939

RR_tauRNAonly_C <- 0.709

RR_eta_C <- 0.62

RR_tauRNA_P <- 0.939


RR_tauRNAonly_P <- 0.709

RR_eta_P <- 0.89

RRlst <- list("C" = list("tau_RNA" = RR_tauRNA_C, 
                         "tau_poct" = RR_tauRNAonly_C,
                         "eta" = RR_eta_C),
              "P" = list("tau_RNA" = RR_tauRNA_P, 
                         "tau_poct" = RR_tauRNAonly_P,
                         "eta" = RR_eta_P))

# setting coverage of national program then calibrate to the coverage of national program treatments

Ccal_C <- 0.01208

Ccal_P <- 0.05424902


Ccal <- list("C" = Ccal_C, "P" = Ccal_P)
# set up parameters for scenario 
dfList_NP <- dfList 
POC_AU$cabY
S_Yint <- 2022
S_Yend <-2024
# function to estimate the parameters of national program 
Param_cal <- function(pj, dlist, index ,S_Yint, S_Yend, r_Yend,ORlist, Ccal ){ 
  #proj: project name 
  # dlist: list of cascade parameters 
  # index: the parameter to calculate 
  # SYint: starting year of scenario 
  # S_Yend: end year of scenario 
  # Orlist: the oddration of the parameters in scenario (save in a list)
  # Ccal: the list of national program coverage 
  # convert to monthly probability then convert back to yearly probability 
  SYpoint_int <- (S_Yint - pj$cabY)/pj$timestep + 1
  SYpoint_end <- (S_Yend - pj$cabY)/POC_AU$timestep
  
  SY_leng <- SYpoint_end - SYpoint_int + 1 
  # remain the same level to end of which year
  rYpoint_end <- (r_Yend - pj$cabY)/POC_AU$timestep
  
  rY_leng <- rYpoint_end - SYpoint_end
  
  intVal <- dlist[[index]][, 3, (SYpoint_int - 1)]
  
  endVal <- c() 
  endVal[1] <- 1-(1- intVal[1])^pj$timestep + (1-(1- (ORlist[["C"]][[index]])*Ccal[["C"]])^pj$timestep)
  endVal[2] <- 1-(1- intVal[2])^pj$timestep + (1-(1- (ORlist[["C"]][[index]])*Ccal[["C"]])^pj$timestep)
  endVal[3] <- 1-(1- intVal[3])^pj$timestep + (1-(1- (ORlist[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
  endVal[4] <- 1-(1- intVal[4])^pj$timestep + (1-(1- (ORlist[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
  endVal[5] <- 1-(1- intVal[5])^pj$timestep + (1-(1- (ORlist[["P"]][[index]])*Ccal[["P"]])^pj$timestep)
  
  endVal[1] <- 1- (1 - endVal[1])^(1/pj$timestep)
  endVal[2] <- 1- (1 - endVal[2])^(1/pj$timestep)
  endVal[3] <- 1- (1 - endVal[3])^(1/pj$timestep)
  endVal[4] <- 1- (1 - endVal[4])^(1/pj$timestep)
  endVal[5] <- 1- (1 - endVal[5])^(1/pj$timestep)


  if(index== "tau_poct"){ 
    avcov <- dlist[["tau_ab"]][, 3, SYpoint_int]
    endVal[1] <- 
      ifelse((endVal[1] + avcov[1])>= 1, (1 - avcov[1]), endVal[1])
    endVal[2] <- 
      ifelse((endVal[2] + avcov[2])>= 1, (1 - avcov[2]), endVal[2])
    endVal[3] <- 
      ifelse((endVal[3] + avcov[3])>= 1, (1 - avcov[3]), endVal[3])
    endVal[4] <- 
      ifelse((endVal[4] + avcov[4])>= 1, (1 - avcov[4]), endVal[4])
    endVal[5] <- 
      ifelse((endVal[5] + avcov[5])>= 1, (1 - avcov[5]), endVal[5])
  }
  
  
  
  for ( i in 2:dim(dlist[[index]])[[2]]){
    dlist[[index]][1, i, c((SYpoint_int - 1):pj$npts)] <- 
      c(seq(as.numeric(intVal[1]), as.numeric(endVal[1]), length = (SY_leng + 1)), 
        rep(as.numeric(endVal[1]), length = (rY_leng)),
        rep(intVal[1], pj$npts - rYpoint_end))
    
    dlist[[index]][2, i, c((SYpoint_int - 1):pj$npts)] <- 
      c(seq(as.numeric(intVal[2]), as.numeric(endVal[2]), length = (SY_leng + 1)), 
        rep(as.numeric(endVal[2]), length = (rY_leng)),
        rep(intVal[2], pj$npts - rYpoint_end))
    
    dlist[[index]][3, i, c((SYpoint_int - 1):pj$npts)] <- 
      c(seq(as.numeric(intVal[3]), as.numeric(endVal[3]), length = (SY_leng + 1)), 
        rep(as.numeric(endVal[3]), length = (rY_leng)),
        rep(intVal[3], pj$npts - rYpoint_end))
    
    dlist[[index]][4, i, c((SYpoint_int - 1):pj$npts)] <- 
      c(seq(as.numeric(intVal[4]), as.numeric(endVal[4]), length = (SY_leng + 1)), 
        rep(as.numeric(endVal[4]), length = (rY_leng)),
        rep(intVal[4], pj$npts - rYpoint_end))
    
    dlist[[index]][5, i, c((SYpoint_int - 1):pj$npts)] <- 
      c(seq(as.numeric(intVal[5]), as.numeric(endVal[5]), length = (SY_leng + 1)), 
        rep(as.numeric(endVal[4]), length = (rY_leng)),
        rep(intVal[4], pj$npts - rYpoint_end))
    
    
  }
  
  return(dlist[[index]])
}

param_var <- c("tau_RNA", "tau_poct", "eta")
# Scenario: current achievement 
dfList_NP <- dfList
for(i in param_var){ 
  dfList_NP[[i]] <- Param_cal(pj = POC_AU, dlist = dfList, index = i, 
                              S_Yint = 2022, S_Yend = 2024, r_Yend = 2024,ORlist = RRlst, 
                              Ccal = Ccal)
}

# Scenario: current achievement + 2024 prediction 
# odds of coverage: 1.08 (total test_predicited/ number test in 2023)
Ccal2024 <- lapply(Ccal,function(x) x*1.08)
dfList_NP2024 <- dfList_NP
for(i in param_var){ 
  dfList_NP2024[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                              S_Yint = 2024, S_Yend = 2026, r_Yend = 2026,ORlist = RRlst, 
                              Ccal = Ccal2024)
}

save(dfList_NP2024, Ccal2024 ,
     file = file.path(OutputFolder ,
                      paste0(project_name,"NP_2024" ,".rda")))

# scenarioA: remaining at level of 2024 to end of 2027 
dfList_NPscale_A <- dfList_NP
for(i in param_var){ 
  dfList_NPscale_A[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                                  S_Yint = 2024, S_Yend = 2026, r_Yend = 2028,ORlist = RRlst, 
                                  Ccal = Ccal2024)
}
# scenarioB: further increase coverage in 2025 
# 2025 coverage: cov_2024*1.5
dfList_NPscale_B <- dfList_NP




# extract the additional probability contributed by national program 
extra_cascade <- function(dList, dList_sc){ 
  param_sc <- lapply(names(dList), function(x) a <- dList_sc[[x]] - dList[[x]])
  names(param_sc) <- names(dList)
  
  return(param_sc)
  }

param_NP <- extra_cascade(dfList, dfList_NP)
param_NP$tau_poct[, , 85]
# import cost data
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
tic <- proc.time()

endY <- 100

param_dfList <- lapply(dfList, function(x) x*0)
names(param_dfList) <- names(dfList)


Sce_sq <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                 disease_progress,pop_array,
                 dfList, param_cascade_sc = param_dfList , fib = fib, 
                 modelrun="UN", proj = "POC_AU", end_Y = endY, 
                 cost = costdfList, costflow = costflow, 
                 costflow_Neg = costflow_Neg)

Sce_np <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                 disease_progress,pop_array,
                 dfList, param_cascade_sc = param_sc , fib = fib, 
                 modelrun="UN", proj = "POC_AU", end_Y = endY, 
                 cost = costdfList, costflow = costflow, 
                 costflow_Neg = costflow_Neg)





toc <- proc.time() - tic 
toc

save(Sce_sq,
     file = file.path(OutputFolder ,
                      paste0(project_name,"Sce_sq" ,".rda")))

save(Sce_np,
     file = file.path(OutputFolder ,
                      paste0(project_name,"Sce_np" ,".rda")))


