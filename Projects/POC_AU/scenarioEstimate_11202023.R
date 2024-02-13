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
Param_cal <- function(pj, dlist, index ,S_Yint, S_Yend, ORlist, Ccal ){ 
  #proj: project name 
  # dlist: list of cascade parameters 
  # index: the parameter to calculate 
  # SYint: starting year of scenario 
  # S_Yend: end year of scenario 
  # Orlist: the oddration of the parameters in scenario (save in a list)
  # Ccal: the list of national program coverage 
  SYpoint_int <- (S_Yint - pj$cabY)/pj$timestep + 1
  SYpoint_end <- (S_Yend - pj$cabY)/POC_AU$timestep + 1
  SY_leng <- SYpoint_end - SYpoint_int
  
  intVal <- dlist[[index]][, 3, SYpoint_int] 
  
  endVal <- c() 
  endVal[1] <- intVal[1] + (ORlist[["C"]][[index]])*Ccal[["C"]]
  endVal[2] <- intVal[2] + (ORlist[["C"]][[index]])*Ccal[["C"]]
  endVal[3] <- intVal[3] + (ORlist[["P"]][[index]])*Ccal[["P"]]
  endVal[4] <- intVal[4] + (ORlist[["P"]][[index]])*Ccal[["P"]] 
  endVal[5] <- intVal[5] 
  
  endVal[1] <- ifelse(endVal[1]>= 1, 1, endVal[1])
  endVal[2] <- ifelse(endVal[2]>= 1, 1, endVal[2])
  endVal[3] <- ifelse(endVal[3]>= 1, 1, endVal[3])
  endVal[4] <- ifelse(endVal[4]>= 1, 1, endVal[4])
  
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
      intVal[5]
  }
  
  
  
  for ( i in 2:dim(dlist[[index]])[[2]]){
    dlist[[index]][1, i, c(SYpoint_int:pj$npts)] <- 
      c(seq(as.numeric(intVal[1]), as.numeric(endVal[1]), length = (SY_leng)), 
        rep(endVal[1], pj$npts - SYpoint_end +1))
    
    dlist[[index]][2, i, c(SYpoint_int:pj$npts)] <- 
      c(seq(as.numeric(intVal[2]), as.numeric(endVal[2]), length = (SY_leng)), 
        rep(endVal[2], pj$npts - SYpoint_end +1))
    
    dlist[[index]][3, i, c(SYpoint_int:pj$npts)] <- 
      c(seq(as.numeric(intVal[3]), as.numeric(endVal[3]), length = (SY_leng)), 
        rep(endVal[3], pj$npts - SYpoint_end +1))
    
    dlist[[index]][4, i, c(SYpoint_int:pj$npts)] <- 
      c(seq(as.numeric(intVal[4]), as.numeric(endVal[4]), length = (SY_leng)), 
        rep(endVal[4], pj$npts - SYpoint_end +1))
    
    
  }
  
  return(dlist[[index]])
}

param_var <- c("tau_RNA", "tau_poct", "eta")

for(i in param_var){ 
  dfList_NP[[i]] <- Param_cal(pj = POC_AU, dlist = dfList_NP, index = i, 
                              S_Yint = 2022, S_Yend = 2024, ORlist = RRlst, 
                              Ccal = Ccal)
}

save(dfList_NP, Ccal,
     file = file.path(OutputFolder ,
                      paste0(project_name,"S_NP_test" ,".rda")))

param_sc <- lapply(names(dfList), function(x) a <- dfList_NP[[x]] - dfList[[x]])
names(param_sc) <- names(dfList)

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


