# this script is for finding the best fitted taegeted coverage of PrEP and HIV patients on treatment 
# once find the targeted coverage of PrEP and HIV diagnosis rate, producing the parameter sets for these scenarios 

rm(list = ls()) 
# library
library(readxl)
library(dplyr)
library(here)
library(data.table)
library(doMC)
library(ggplot2)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library(doParallel)
registerDoMC(cores = 2) 

# set up dircetory 
codep <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/TWHCV-model"
epidatapath <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Taiwan-MSM-HCV-model"
#file path of "TWHCV-model" project
## !! using the github simplfied HCV testing v1.0.1 release on 2023.04.17
Rcode <- file.path(codep, '03. Code')
scenariopath <- file.path(epidatapath, '01. DATA/model input/Scenarios')

# save and cost data folder 
dtp <- "/Users/jjwu/Library/CloudStorage/OneDrive-UNSW/05. PhD Project/Simplified HCV testing model_/Projects/02.CEA_TWHCVMSM"
DataFolder <- file.path(dtp, "01. DATA/model input")
outputdt <- file.path(dtp, "02. Output/RDA")
# Rda file path 
# load the .rda file of base estimate 


load(file.path(epidatapath , paste0("HCVModel", "param",".rda")))

load(file.path(epidatapath , paste0("HCVModel",".rda")))



# source 
source(file.path(Rcode, "Functions/HCV_model.R"))

source(file.path(Rcode, "Functions/plotFunctions.R"))

source(file.path(Rcode, "Functions/plotOptions.R"))

source(file.path(Rcode, "Functions/Scenarios_wrangling.R"))



load(file.path(epidatapath , paste0("HCVModel",".rda")))

load(file.path(epidatapath, paste0("HCVModel", "cali",".rda")))

load(file.path(outputdt, "cost.rda"))

load(file.path(epidatapath, paste0("HCVModel", "param",".rda")))

# output file path 
# dir.create("02. Output") # create subdircetory 
# dir.create("02. Output/RDA") 


# source 
source(file.path(Rcode, "Functions/HCV_model.R"))

source(file.path(Rcode, "Functions/plotFunctions.R"))

source(file.path(Rcode, "Functions/plotOptions.R"))

source(file.path(Rcode, "Functions/Scenarios_wrangling.R"))

# inititaing year 
calibrateY <- 2004

scenario_Y <- 2022

uptake_Y <- 2024 

PrEPcov_Y <- 2030 

# uptake PrEP rate 

uPrEP <- 15
uHIVd <- 2.05
# finding time point of start of 2022
scenario_pt <- length(seq(calibrateY, scenario_Y, by = HCV$timestep))

uptake_pt <- length(seq(calibrateY, uptake_Y, by = HCV$timestep))

PrEPcov_pt <- length(seq(calibrateY, PrEPcov_Y, by = HCV$timestep))

PrEPcov_pt_off <- length(seq(calibrateY, PrEPcov_Y + 1, by = HCV$timestep))

# uptake PrEP
pop_arrayPrEPHIV <- list(pop_array, pop_array,pop_array)

pop_arrayPrEPHIV[[1]][1,2,c(scenario_pt:PrEPcov_pt)] <-  
  
  seq(as.numeric(pop_array[1, 2, scenario_pt]), 
      
      uPrEP*as.numeric(pop_array[1, 2, scenario_pt]), 
      
      length = (PrEPcov_pt - scenario_pt + 1))

pop_arrayPrEPHIV[[1]][1,2,c(PrEPcov_pt: dim(pop_arrayPrEPHIV[[1]])[3])] <-  
  
  uPrEP*as.numeric(pop_array[1, 2, scenario_pt])


# uptake HIV diagnosis rate

pop_arrayPrEPHIV[[2]][3,4,c(scenario_pt:PrEPcov_pt)] <-  
  
  seq(as.numeric(pop_array[3, 4, scenario_pt]), 
      
      uHIVd*as.numeric(pop_array[3, 4, scenario_pt]), 
      
      length = (PrEPcov_pt - scenario_pt + 1))

pop_arrayPrEPHIV[[2]][3,4,c(PrEPcov_pt: dim(pop_arrayPrEPHIV[[2]])[3])] <-  
  
  uHIVd*as.numeric(pop_array[3, 4, scenario_pt])


# uptake HIV diagnosis and PrEP 

pop_arrayPrEPHIV[[3]][3,4,] <-  pop_arrayPrEPHIV[[2]][3,4,]

pop_arrayPrEPHIV[[3]][1,2,] <- pop_arrayPrEPHIV[[1]][1,2,]

names(pop_arrayPrEPHIV) <- c("PrEP", "HIVD", "PrEPnHIVD")

# we specify the number of cores/workers we want to use
tic <- proc.time()

# Run model on best estimates

Sq <- list()
for(i in names(pop_arrayPrEPHIV)){
  Sq[[i]] <- HCVMSM(HCV,best_estimates, best_est_pop,
                        disease_progress,pop_arrayPrEPHIV[[i]], dfList, fib,
                        end_Y = 30, modelrun="UN", cost = cost, 
                        costflow = cost$flow)
  }


# HIV diag %
options(warn=-1)

HIVdiag <- read.csv(file.path(paste0(rdapath, "/01. DATA", sep="/"), 
                              "HIVdiag.csv"), header = TRUE)%>%
  as.data.frame()%>%
  mutate(time = year - 2003, realPop = HIV.diagnosis.rate*100,
         low = lower*100,
         up = upper*100)

# subtract the number of diagnosed in the model
HIVInfectedD <- list()

HIVDiagRate <- list()

HIVInfected <- list() 

HIVInf <- list()

HIVN <- list()

HIVPrEPRate <- list()

HIVPrEP <- list() 

HIVNeg <- list()

for(i in seq_along(Sq)){ 
  
  HIVInfected[[i]] <- popResults_MidYearBase(HCV, Sq[[i]],
                                         Population = c("HIV+", "HIV+d"),
                                         Disease_prog = NULL, 
                                         Cascade = NULL, param = NULL, 
                                         endYear = 30) 
  
  HIVInf[[i]] <- HIVInfected[[i]]%>%ungroup()%>% dplyr::group_by(year)%>%
    summarise(best = sum(best))
  
  HIVInfectedD[[i]] <- popResults_MidYearBase(HCV, Sq[[i]],Population = c("HIV+d"),
                                          Disease_prog = NULL, 
                                          Cascade = NULL, param = NULL,
                                          endYear = 30)
  HIVDiagRate[[i]] <- cbind(year = seq(HCV$startYear , 30-1 ,1),
                            as.data.frame(HIVInfectedD[[i]][, -c(1,2)]/ HIVInf[[i]][ ,-1])*100)%>%
    tibble::as_tibble()%>%mutate(year = year + HCV$cabY - 1)
  
  HIVN[[i]] <- popResults_MidYearBase(HCV, Sq[[i]],
                                  Population = c("HIV-", "HIV-PrEP"),
                                  Disease_prog = NULL, 
                                  Cascade = NULL, param = NULL, 
                                  endYear = 30) 
  
  HIVNeg[[i]] <- HIVN[[i]]%>%ungroup()%>% dplyr::group_by(year)%>%
    summarise(best = sum(best))
  
  HIVPrEP[[i]] <- popResults_MidYearBase(HCV, Sq[[i]],Population = c("HIV-PrEP"),
                                     Disease_prog = NULL, 
                                     Cascade = NULL, param = NULL,
                                     endYear = 30)
  HIVPrEPRate[[i]] <- cbind(year = seq(HCV$startYear , 30-1 ,1),
                            as.data.frame(HIVPrEP[[i]][, -c(1,2)]/HIVNeg[[i]][ ,-1])*100)%>%
    tibble::as_tibble()%>%mutate(year = year + HCV$cabY - 1)
  
 }
Diag <- list()
Pr <- list()
for(i in seq_along(Sq)){ 
  
  Diag[[i]] <- HIVDiagRate[[i]]%>%filter(year == 2030)
  
  Pr[[i]] <- HIVPrEPRate[[i]]%>%filter(year == 2030)
  
  }

# uncertainty range for these best fit 
## trim the timepoints 
trimpt <- 1000
#### increasing PrEP uptake & HIV diagnosis from 2022-2030 ####

# double the PrEP rate 
Param_pop_arrayPrEPHIV <- list()

for(i in names(pop_arrayPrEPHIV)){ 
  
  pop_arrayPrEPHIV[[i]] <- pop_arrayPrEPHIV[[i]][ , ,c(1:trimpt)]
  
  Param_pop_arrayPrEPHIV[[i]] <- lapply(randomParams_array, function(x) { 
    
    a <- x*pop_arrayPrEPHIV[[i]]
  
    })

  
  }

costTrim <- list()
for(i in names(cost)){ 
  
  costTrim[[i]] <- cost[[i]][, , c(1:trimpt)]
  }


Param_costList <- lapply(names(costTrim), function(m) lapply(
    randomParams_array, function(x){ x * costTrim[[m]] }))

Param_costList  <- purrr::transpose(Param_costList)

for(i in 1:length(Param_costList)){
  names(Param_costList[[i]]) <- names(costTrim)
  
}

save(Param_costList, costTrim,
     file = file.path(outputdt, "cost.rda"))

save(pop_arrayPrEPHIV, Param_pop_arrayPrEPHIV,
     file = file.path(outputdt, "mainS_PrEPHIV.rda"))


#### turn tranmission off at the end of 2030 
best_estimatesOff <- best_estimates

best_estimatesOff[c(PrEPcov_pt_off: nrow(best_estimatesOff )), "beta1"] <- 0

best_estimatesOff[c(PrEPcov_pt_off: nrow(best_estimatesOff )), "beta2"] <- 0

best_estimatesOff[c(PrEPcov_pt_off: nrow(best_estimatesOff )), "beta3"] <- 0

best_estimatesOff[c(PrEPcov_pt_off: nrow(best_estimatesOff )), "beta4"] <- 0

Param_estimatesOff <- Param_estimates

for(set in 1:HCV$numberSamples){ 
  
  Param_estimatesOff[[set]][c(PrEPcov_pt_off: nrow(best_estimatesOff )), "beta1"] <- 0
  
  Param_estimatesOff[[set]][c(PrEPcov_pt_off: nrow(best_estimatesOff )), "beta2"] <- 0
  
  Param_estimatesOff[[set]][c(PrEPcov_pt_off: nrow(best_estimatesOff )), "beta3"] <- 0
  
  Param_estimatesOff[[set]][c(PrEPcov_pt_off: nrow(best_estimatesOff )), "beta4"] <- 0
  
  Param_estimatesOff[[set]] <- setnafill(Param_estimatesOff[[set]], type = "locf") # fill na value with last row in each column
  
  }





save(best_estimatesOff, Param_estimatesOff, file = file.path(outputdt, "toft.rda"))
