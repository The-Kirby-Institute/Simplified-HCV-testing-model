# This script including the time-varying parameters adjust
# which needs to be revised to idea format as 1) import the csv file for 
#                                                the parameter changing value at specify year
#                                             2) using exponential to revise "best_estimate"
# Reorganzed also required in this script
# To seperate the summary results 
rm(list = ls())
## manually calibration 
library(ggplot2)
library(ggrepel)
library(directlabels)
library(data.table)
library(gridExtra)
library(grid)
basePath <- getwd()
projectFolder <- file.path(basePath)
project_name <- "HCVModel"
DataFolder <- file.path(basePath, "01. DATA/model input" )
Rcode <- file.path(basePath, "03. Code")
urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- TRUE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

source(file.path(Rcode, "/Functions/HCV_model.R"))
source(file.path(Rcode, "/Functions/plotFunctions.R"))

load(file.path(projectFolder, paste0(project_name, ".rda")))
load(file.path(projectFolder, paste0(project_name, "cali", ".rda"))) 

# % increasing for stages in cascade

HCVab <- 0.25
HCVab_hivd <- 0.8
HCVRNA_hivd <- 0.8
TreatInit <- 0.26989
TreatInit_hivd <- 0.678
retreat <- 0.1
SVR <- 0.919

## entry
#### manually change the entry number  
 MSMentry <- read.csv(file.path(paste0(projectFolder, "/01. DATA", sep="/"),
                               "MSM_entry.csv"), header = TRUE)%>%
  as.data.frame()

 
 remainY <- last(MSMentry$year) + 1 - HCV$cabY + 1 
 
 MSMentry[c(remainY:HCV$endYear), 1] <- seq(last(MSMentry$year) + 1, 
                                            last(MSMentry$year) + 1 + 
                                              (HCV$endYear - remainY), 1)

 MSMentry2021 <- MSMentry[c(1:remainY - 1 ), 2]

 MSMafrer2021 <- rep(MSMentry[(remainY - 1), 2], HCV$endYear - remainY + 1)

 MSM_Entry <- MSMentry%>%
  mutate(number = c(MSMentry2021, MSMafrer2021))%>%
  mutate(number = number * best_estimates$MSM_pro[1])

 Ent <- as.data.frame(matrix(0, ncol = 1,
                            nrow = HCV$npts)) 


 for (year in 1:(HCV$nyears - 1)){
  
  indices <- ((year - 1) * 1/HCV$timestep + 1): (year * 1/HCV$timestep + 1)
  
  yearvalue <- seq(MSM_Entry$number[year], MSM_Entry$number[year + 1], 
                   length = (1 + 1/HCV$timestep))
  
  Ent[indices, 1] <-yearvalue 
  
 }


 best_estimates$entry <- Ent$V1  

#### manually change beta value  

#### 2004 - 2008 fixed as initial beta value 
#### increasing beta from 2008 -2018 
#### 2018 - fixed 


 calibration_Y <- 2004 
 varying_Y <- c(2014, 2014, 2014, 2014)
 
 varyingYpoint <- (varying_Y - 2004)/HCV$timestep + 1
 #varyingYpoint2 <- (2018 - varying_Y)/HCV$timestep+1
 #varyingYpoint3 <- (2025 - 2018)/HCV$timestep+1
 betaVaring <- best_estimates%>%select(beta1, beta2, beta3, beta4)%>%
   slice(1)
 
 # changing beta value 
 betaVaring[2,] <- c(0.08, 0.08, 1.1, 1)
 
 #betaVaring[3,] <- c(0.08, 0.08, 1.2, 1.2)
 # betaVaring[2,] <- c(0.2, 0.15, 2.2, 1.2)
 # different scenario: beta keep increasing to 2030
 # slope 
 # slp <- (betaVaring[2,] - betaVaring[1,])/varyingYpoint 
 
 # var_Y <- 2025
 # var_Ypoint <- (var_Y-calibration_Y)/HCV$timestep+1
 # betaVaring[2, ] <- betaVaring[1, ] + slp*var_Ypoint
 
 # append the new beta at specific time point 
 # var_Ypoint <- (var_Y-calibration_Y)/HCV$timestep+1
 
 # var_Ypoint_PrEP <- (var_Y-2017)/HCV$timestep+1
 # PrEPbeta <- betaVaring[1,2 ] + slp[1]*var_Ypoint_PrEP
 # betaVaring[2, 2] <- PrEPbeta
 yearvalue <- as.data.frame(matrix(0, ncol = ncol(betaVaring),
                                   nrow = HCV$npts + 1))  
 yearvalue[, 1] <- c(seq(betaVaring[1, 1], betaVaring[2, 1], 
                         length = varyingYpoint[1]), 
                     rep(betaVaring[2, 1], (HCV$npts - varyingYpoint[1]) + 1))
 
 
 # start from 2017
 
 # 2004-2017: 0
 
 PrEP_Y <- 2017
 var_before <- (PrEP_Y - calibration_Y)/HCV$timestep - 1
 # value for PrEP users after 2017 
 
 yearvalue[, 2] <- c(rep(0, var_before), rep(betaVaring[2, 2], 
                                            HCV$npts - var_before + 1))
 yearvalue[, 3] <- c(seq(betaVaring[1, 3], betaVaring[2, 3], 
                         length = varyingYpoint[3]), 
                     rep(betaVaring[2, 3], (HCV$npts - varyingYpoint[3] ) + 1))
 
 yearvalue[, 4] <- c(seq(betaVaring[1, 4], betaVaring[2, 4], 
                         length = varyingYpoint[4]), 
                     rep(betaVaring[2, 4], (HCV$npts - varyingYpoint[4] ) + 1))
 
 best_estimates$beta1 <- yearvalue[, 1]
 best_estimates$beta2 <- yearvalue[, 2]
 best_estimates$beta3 <- yearvalue[, 3]
 best_estimates$beta4 <- yearvalue[, 4]
 
 
 # HIV testing rate from HIV+ to HIV + d
 
 varying_Y <- 2020
 varyingYpoint <- (varying_Y-calibration_Y)/HCV$timestep+1
 
 pop_array[3,4, ] <- c(seq(0.28, 0.28,length = varyingYpoint), 
                       rep(0.28, HCV$npts - varyingYpoint))
 

# considering introduced of DAA, parameter change from 2019 to 2021   
# eta: treatment initiated rate change since 2019 regarding to DAA era

varying_Yint <- 2017
varying_Yend <-2021
varyingYpoint_int <- (varying_Yint - calibration_Y)/HCV$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/HCV$timestep + 1


# HCV antibody testing 

intVal <- dfList$tau_ab[1, 3, 1] 

for ( i in 2:dim(dfList$eta)[[2]]) {
  for ( j in 1: (HCV$npops - 1)){
    dfList$tau_ab[j, i, c(varyingYpoint_int:HCV$npts)] <- 
      c(seq(intVal, HCVab, length = varyingYpoint_end - varyingYpoint_int), 
        rep(HCVab, HCV$npts - varyingYpoint_end + 1)) 
    
    
  }
  
}

intVal <- dfList$tau_ab[4, 3, 1] 
for ( i in 2:dim(dfList$eta)[[2]]){
    dfList$tau_ab[4, i, c(varyingYpoint_int:HCV$npts) ] <- 
      c(seq(intVal, HCVab_hivd , length = varyingYpoint_end - varyingYpoint_int), 
        rep(HCVab_hivd , HCV$npts - varyingYpoint_end + 1)) 
    
    
  }

 intVal <- dfList$tau_ag[4, 3, 1] 
 for ( i in 2:dim(dfList$tau_ag)[[2]]){
  dfList$tau_ag[4, i, c(varyingYpoint_int:HCV$npts) ] <- 
    c(seq(intVal, HCVRNA_hivd, length = varyingYpoint_end - varyingYpoint_int), 
      rep(HCVRNA_hivd, HCV$npts - varyingYpoint_end + 1)) 
  
  
 }

# coverage to every time step 

# intial value 
intVal <- dfList$eta[1, 3, 1] 

for ( i in 2:dim(dfList$eta)[[2]]){
  dfList$eta[c(1:3), i, c(varyingYpoint_int:HCV$npts)] <- 
    c(seq(intVal, TreatInit, length = varyingYpoint_end - varyingYpoint_int),
      rep(TreatInit, HCV$npts - varyingYpoint_end + 1))
  
  
    dfList$eta[4, i, c(varyingYpoint_int:HCV$npts)] <- 
      c(seq(intVal, TreatInit_hivd, length = varyingYpoint_end - varyingYpoint_int),
       rep(TreatInit_hivd, HCV$npts - varyingYpoint_end + 1)) 
    
    
  }


# lota: treatment failed
# we assumed it is not changed after introduced of DAA regarding to 
# there were HCV case managers to follow-up patients situation


#　rho:re-treated
intVal <- dfList$rho[1, 4, 1] 

for ( i in 2:dim(dfList$eta)[[2]]){
  for ( j in 1:HCV$npops){
    dfList$rho[j, i, c(varyingYpoint_int:HCV$npts)] <- 
      c(seq(intVal, retreat, length = varyingYpoint_end - varyingYpoint_int), 
        rep(retreat, HCV$npts - varyingYpoint_end + 1)) 
    
    
  }
}

#SVR: cure  
intVal <- dfList$cured[1, 2, 1] 

for ( i in 2:(dim(dfList$cured)[[2]] - 1)){
  for ( j in 1: HCV$npops){
    dfList$cured[j, i, c(varyingYpoint_int:HCV$npts) ] <- 
      c(seq(intVal, SVR, length = varyingYpoint_end - varyingYpoint_int), 
        rep(SVR, HCV$npts - varyingYpoint_end + 1)) 
    
    
  }
}
save(HCV, constants, disease_progress, fib, dfList, pop_array, 
     constantsDf, initialPops, best_estimates, best_initial_pop, 
     file = file.path(projectFolder,
                      paste0(project_name, ".rda")))


