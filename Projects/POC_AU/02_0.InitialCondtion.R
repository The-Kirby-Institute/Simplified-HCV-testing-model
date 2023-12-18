# Running 02_1.CalibrateOptimze.R and revised the input parameters then 
# Run 01.SetupModel to align the input parameters to fit the format that model need 
# Then Run this script to get the initial condition 

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
# load(file.path(OutputFolder, paste0(project_name, "cali.rda")))


urrTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
runSamples <- FALSE
saveAsBase <- TRUE  # if TRUE doesn't append time to results and overwrites
# a base file. Useful for storing main results or
# testing

source(file.path(Rcode, "/Functions/HCV_model.R"))

source(file.path(Rcode, "/Functions/plotFunctions.R")) 
source(file.path(Rcode, "/Functions/check_steady.R"))

pop_transitions <- read.csv(file.path(DataFolder, "population_transitions.csv"), 
                            header = TRUE)

pop_transitions <- pop_transitions[ , -1] # drop the row number

# get first row only 
pop_transitions <- pop_transitions[1,]

xx <- array(unlist(pop_transitions), c(POC_AU$npops, POC_AU$npops, 1))

xx_array <- aperm(xx, c(2, 1, 3))

colnames(xx_array) <- POC_AU$popNames

rownames(xx_array) <- POC_AU$popNames

xx_array["C_PWID","C_fPWID", 1] <- 0.065
xx_array["C_PWID","P_PWID", 1] <- 0.115
xx_array["C_fPWID","C_PWID", 1] <- 0.0015
xx_array["C_fPWID","P_fPWID", 1] <- 0.021
xx_array["P_PWID","C_PWID", 1] <- 0.9
xx_array["P_PWID","P_fPWID", 1] <- 0.5
xx_array["P_fPWID","C_fPWID", 1] <- 0.82
xx_array["P_fPWID","P_PWID", 1] <- 0.08



# extend to same length 
pop_array  <- array(matrix(xx_array[ , ,1]), 
                    c(POC_AU$npops, POC_AU$npops,POC_AU$npts))


best_estimates$beta1 <- 0.3

best_estimates$beta2 <- 0.08
best_estimates$beta3 <- 1.3
best_estimates$beta4 <- 0.35
best_estimates$beta5 <- 0.08
best_estimates$HCVP1 <- 0.48
best_estimates$HCVP2 <- 0.25
best_estimates$HCVP3 <- 0.9
best_estimates$HCVP4 <- 0.3
best_estimates$HCVP5 <- 0.003



# the current best-estimated 

#save(best_estimates,pop_array, 
#     file = file.path(OutputFolder , paste0(project_name,"temp_best" ,".rda")))


tic <- proc.time()

steady <- HCVMSM(POC_AU, best_estimates, best_initial_pop,
                 disease_progress, pop_array,
                 dfList, fib, end_Y = 1000, 
                 modelrun = "steady", proj = "POC_AU")

toc <- proc.time() - tic 


toc



#### extract infected proportion as the infected population allocation #### 
# finding equilibrium

check_steady(model_result = steady, endY = 1000,
             timestep = POC_AU$timestep, 
             Ncomp = POC_AU$ncomponent*POC_AU$npops, 
             Tequilibrium = 1000)


#extract proportion of pops for initial condition 

df_list <- lapply(steady, as.data.frame.table)

allpop <- df_list$allPops%>%mutate(time = rep(seq(1,(1000 - POC_AU$timestep),POC_AU$timestep), 
                                              each=POC_AU$ncomponent*POC_AU$npops),
                                   Frequency=Freq)

popPro_extract <- df_list$allPops%>%
  mutate(time = rep(seq(POC_AU$startYear, (1000 - POC_AU$timestep), POC_AU$timestep), 
                    each=POC_AU$ncomponent * POC_AU$npops),
         Frequency=Freq)%>%
  filter(time == 500)%>%
  mutate(cascade_status = sub("^[^_]*_", "", Var2), 
         dis_prog = sub("\\_.*", "", Var2),
         SI = ifelse(cascade_status%in%c("s", "cured"), "S","I"),
         parameter =Var2)%>%group_by(Var1 ,SI)%>%
  mutate(total = sum(Frequency),
         value = ifelse(Frequency==0, 0, Frequency/total))%>%
  ungroup()%>%group_by(Var1)%>%mutate(pop_prop = ifelse(
    Frequency==0, 0, Frequency/sum(Frequency)))%>%
  ungroup()%>%dplyr::select(Var1,parameter, value, SI)


write.csv(popPro_extract, 
          file.path(DataFolder,"/Estimate_initial_pop.csv")) 

#### number of PWID/former PWID in each population ####

estPops <- popPro_extract


init_pop <- filter(initialPops, parameter == "init_pop")$value

pop_prop <- initialPops%>%filter(parameter%in% c("pop_prop1", "pop_prop2", 
                                                 "pop_prop3", "pop_prop4", 
                                                 "pop_prop5"))%>%
  dplyr::select(value)%>%unlist()%>%as.vector()

popProp <- as.numeric(init_pop)*pop_prop 


# prevalence at initial
#init_prop_S <- initialPops%>%filter(parameter%in% c("s1", "s2", "s3", "s4", "s5"))%>%
#  select(value)%>%unlist()%>%as.vector()
init_prop_I <- c(best_estimates$HCVP1[1], best_estimates$HCVP2[1], 
                 best_estimates$HCVP3[1], best_estimates$HCVP4[1],
                 best_estimates$HCVP5[1])
init_prop_S <- c(1 - init_prop_I)


best_estimates$HCVP1

estPops <- estPops%>%mutate(
  pop_group = rep(c(popProp),dim(estPops)[1]/POC_AU$npops),
  SIprop = ifelse(estPops$SI=="S", 
                  rep(init_prop_S, POC_AU$diseaseprogress_n*POC_AU$npops),
                  rep(init_prop_I, POC_AU$ncomponent*POC_AU$npops - 
                        POC_AU$diseaseprogress_n*POC_AU$npops)),
  est_pop = value*pop_group*SIprop)

best_est_pop <- as.matrix(as.data.frame(matrix(estPops$est_pop, 
                                               ncol = POC_AU$ncomponent,  
                                               nrow = POC_AU$npops)))

colnames(best_est_pop) <- c(POC_AU$component_name)

save(POC_AU,steady, best_est_pop, best_estimates, pop_array,
     file = file.path(OutputFolder ,
                      paste0(project_name,"cali" ,".rda")))


# calibration 
# adding the time-varying 
# timevarying 
# calibration 
# adding the time-varying 
varying_Yint <- 2015
varying_Yend <-2020
calibration_Y <- 2015
varyingYpoint_int <- (varying_Yint - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/POC_AU$timestep + 1


# treatment init 
# coverage to every time step 
#bvalue <- c(0.12, 0.1, 0.95, 0.3, 0.08)
bvalue <- c(best_estimates$beta1[1], best_estimates$beta2[1], 
            best_estimates$beta3[1], best_estimates$beta4[1], 
            best_estimates$beta5[1])

best_estimates$beta1 <- 
  c(seq(best_estimates$beta1[1], bvalue[1], length = (varyingYpoint_end)), 
    rep(bvalue[1], POC_AU$npts - varyingYpoint_end))

best_estimates$beta2 <- 
  c(seq(best_estimates$beta2[1], bvalue[2], length = (varyingYpoint_end)), 
    rep(bvalue[2], POC_AU$npts - varyingYpoint_end))

best_estimates$beta3 <- 
  c(seq(best_estimates$beta3[1], bvalue[3], length = (varyingYpoint_end)), 
    rep(bvalue[3], POC_AU$npts - varyingYpoint_end))

best_estimates$beta4 <- 
  c(seq(best_estimates$beta4[1], bvalue[4], length = (varyingYpoint_end)), 
    rep(bvalue[4], POC_AU$npts - varyingYpoint_end))

TreatInit <- c(0.48,
               0.24,
               0.999999,
               0.999999,
               0.999999)


varying_Yint <- 2015
varying_Yfir <- 2016
varying_Ymid <- 2018
varying_Yend <- 2020
varying_Yend2022 <- 2022
calibration_Y <- 2015
varyingYpoint_int <- (varying_Yint - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_fir <- (varying_Yfir - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_mid <- (varying_Ymid - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end2022 <- (varying_Yend2022 - calibration_Y)/POC_AU$timestep + 1
# intial value 
intVal <- dfList$eta[, 3, 1] 

intVal_rho <- dfList$rho[, 3, 1]

Retreat <- c(intVal_rho[1], intVal_rho[2], 0.999999, 0.999999, 0.999999)

endY <- 10

for ( i in 2:dim(dfList$eta)[[2]]){
  dfList$eta[1, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[1]),as.numeric(TreatInit[1]) , length = (varyingYpoint_fir)),
      seq(as.numeric(TreatInit[1]), as.numeric(TreatInit[1]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[1]), 0.25, length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(0.25, POC_AU$npts - varyingYpoint_end))
  
  dfList$eta[2, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[2]), as.numeric(TreatInit[2]), length = (varyingYpoint_fir)),
      seq(as.numeric(TreatInit[2]), as.numeric(TreatInit[2]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[2]), 0.09, length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(0.09, POC_AU$npts - varyingYpoint_end))
}
varying_Yint <- 2015
varying_Yfir <- 2019
varying_Ymid <- 2020
varying_Yend <- 2021
calibration_Y <- 2015
varyingYpoint_int <- (varying_Yint - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_fir <- (varying_Yfir - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_mid <- (varying_Ymid - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/POC_AU$timestep + 1

for ( i in 2:dim(dfList$eta)[[2]]){  
  dfList$eta[3, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[3]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, as.numeric(TreatInit[3]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[3]), as.numeric(TreatInit[3]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(TreatInit[3]), POC_AU$npts - varyingYpoint_end))
  
  dfList$eta[4, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[4]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, as.numeric(TreatInit[4]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[4]), as.numeric(TreatInit[4]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(TreatInit[4]), POC_AU$npts - varyingYpoint_end))
  
  
  
  dfList$eta[5, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[5]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, as.numeric(TreatInit[5]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[5]), as.numeric(TreatInit[5]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(TreatInit[5]), POC_AU$npts - varyingYpoint_end))
  
  
  dfList$rho[3, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_rho[3]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, as.numeric(Retreat[3]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(Retreat[3]), as.numeric(Retreat[3]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(Retreat[3]), POC_AU$npts - varyingYpoint_end))
  
  dfList$rho[4, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_rho[4]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, as.numeric(Retreat[4]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(Retreat[4]), as.numeric(Retreat[4]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(Retreat[4]), POC_AU$npts - varyingYpoint_end))
  
  
  
  dfList$rho[5, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_rho[5]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, as.numeric(Retreat[5]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(Retreat[5]), as.numeric(Retreat[5]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(Retreat[5]), POC_AU$npts - varyingYpoint_end))
  
}




tic <- proc.time()
endY <- 100
calibrateInit <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                        disease_progress,pop_array,
                        dfList, fib, 
                        modelrun="UN", proj = "POC_AU", end_Y = endY)


toc <- proc.time() - tic 
toc



save(calibrateInit, bvalue, TreatInit,best_estimates, dfList,
     file = file.path(OutputFolder ,
                      paste0(project_name,"cali_timev" ,".rda")))


