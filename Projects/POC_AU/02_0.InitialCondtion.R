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


best_estimates$beta1 <- 0.21
best_estimates$beta2 <- 0.038
best_estimates$beta3 <- 3
best_estimates$beta4 <- 0.3
best_estimates$beta5 <- 0.08
best_estimates$HCVP1 <- 0.5
best_estimates$HCVP2 <- 0.2
best_estimates$HCVP3 <- 0.5
best_estimates$HCVP4 <- 0.3
best_estimates$HCVP5 <- 0.003


# the current best-estimated 

#save(best_estimates,pop_array, 
#     file = file.path(OutputFolder , paste0(project_name,"temp_best" ,".rda")))
dfList_sc <- lapply(dfList, function(x) x*0)
names(dfList_sc) <- names(dfList)
fc <- matrix(0, ncol = dim(dfList$eta)[3], nrow = POC_AU$npops)
tic <- proc.time()

steady <- HCVMSM(POC_AU, best_estimates, best_initial_pop,
                 disease_progress, pop_array,
                 dfList, param_cascade_sc = dfList_sc , fib, end_Y = 1000, 
                 modelrun = "steady", proj = "POC_AU", fc_sc = fc,
                 fp = NULL)

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

TreatInit <- c(0.40,
               0.45,
               0.999999,
               0.999999,
               0.999999)

dfList$eta[, , 1]
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
      seq(as.numeric(TreatInit[1]), as.numeric(TreatInit[1]), length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(as.numeric(TreatInit[1]), POC_AU$npts - varyingYpoint_end))
  
  dfList$eta[2, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[2]), as.numeric(TreatInit[2]), length = (varyingYpoint_fir)),
      seq(as.numeric(TreatInit[2]), as.numeric(TreatInit[2]), length = (varyingYpoint_mid - varyingYpoint_fir)),
      seq(as.numeric(TreatInit[2]), 0.25, length = (varyingYpoint_end - varyingYpoint_mid)),
      rep(0.25, POC_AU$npts - varyingYpoint_end))
  
  
  
}
# ab and rna testing in community increase 

varying_Yint <- 2015
varying_Yfir <- 2019
varying_Ymid <- 2020
varying_Ymid2 <- 2021
varying_Yend <- 2022
calibration_Y <- 2015
varyingYpoint_int <- (varying_Yint - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_fir <- (varying_Yfir - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_mid <- (varying_Ymid - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_mid2 <- (varying_Ymid2 - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/POC_AU$timestep + 1
# intial value 
intVal_ab <- dfList$tau_ab[c(1:2), 3, 1] 
intVal_RNA <- dfList$tau_RNA[c(1:2), 3, 1] 

dfList$tau_ab[, , 61]
endY <- 10

for ( i in 2:dim(dfList$eta)[[2]]){
  dfList$tau_ab[1, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_ab[1]),as.numeric(intVal_ab[1]) , length = (varyingYpoint_fir - 2)),
      seq(as.numeric(intVal_ab[1]), as.numeric(intVal_ab[1]), length =(varyingYpoint_mid - varyingYpoint_fir - 1)),
      seq(as.numeric(intVal_ab[1]), as.numeric(intVal_ab[1]), length = varyingYpoint_mid2 - varyingYpoint_mid + 1),
      seq(as.numeric(intVal_ab[1]), as.numeric(intVal_ab[1])*1.18, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(as.numeric(intVal_ab[1])*1.18), POC_AU$npts - varyingYpoint_end + 1))

  dfList$tau_ab[2, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_ab[2]),as.numeric(intVal_ab[2]) , length = (varyingYpoint_fir - 2)),
      seq(as.numeric(intVal_ab[2]), as.numeric(intVal_ab[2]), length =(varyingYpoint_mid - varyingYpoint_fir - 1)),
      seq(as.numeric(intVal_ab[2]), as.numeric(intVal_ab[2]), length = varyingYpoint_mid2 - varyingYpoint_mid + 1),
      seq(as.numeric(intVal_ab[2]), as.numeric(intVal_ab[2])*1.18, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(as.numeric(intVal_ab[1])*1.18*1.05), POC_AU$npts - varyingYpoint_end + 1))

  dfList$tau_RNA[1, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_RNA[1]),as.numeric(intVal_RNA[1]) , length = (varyingYpoint_mid2 - 2)),
      seq(as.numeric(intVal_RNA[1]), as.numeric(intVal_RNA[1])*1.08, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(intVal_RNA[1])*1.08, POC_AU$npts - varyingYpoint_end + 1))
  
  dfList$tau_RNA[2, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_RNA[2]),as.numeric(intVal_RNA[2]) , length = (varyingYpoint_mid2 - 2)),
      seq(as.numeric(intVal_RNA[2]), as.numeric(intVal_RNA[2])*1.08, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(intVal_RNA[2])*1.08, POC_AU$npts - varyingYpoint_end + 1))

  
}

# ab and rna testing in prison increase 

varying_Yint <- 2015
varying_Yfir <- 2019
varying_Ymid <- 2020
varying_Ymid2 <- 2022
varying_Yend <- 2023
calibration_Y <- 2015
varyingYpoint_int <- (varying_Yint - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_fir <- (varying_Yfir - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_mid <- (varying_Ymid - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_mid2 <- (varying_Ymid2 - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/POC_AU$timestep + 1
# intial value 
intVal_ab <- dfList$tau_ab[c(3:5), 3, 1] 
intVal_RNA <- dfList$tau_RNA[c(3:5), 3, 1] 

dfList$tau_ab[, , 61]
endY <- 10

for ( i in 2:dim(dfList$eta)[[2]]){
  dfList$tau_ab[3, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_ab[1]),as.numeric(intVal_ab[1]) , length = (varyingYpoint_fir - 2)),
      seq(as.numeric(intVal_ab[1]), as.numeric(intVal_ab[1])*1.18, length =(varyingYpoint_mid - varyingYpoint_fir + 1)),
      rep(as.numeric(intVal_ab[1])*1.18, as.numeric(intVal_ab[1])*1.18, length =(varyingYpoint_mid2 - varyingYpoint_mid -1 )),
      seq(as.numeric(intVal_ab[1])*1.18, as.numeric(intVal_ab[1])*1.18*1.05, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(intVal_ab[1])*1.18*1.05, POC_AU$npts - varyingYpoint_end + 1))
  
  dfList$tau_ab[4, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_ab[2]),as.numeric(intVal_ab[2]) , length = (varyingYpoint_fir - 2)),
      seq(as.numeric(intVal_ab[2]), as.numeric(intVal_ab[2])*1.18, length =(varyingYpoint_mid - varyingYpoint_fir + 1)),
      rep(as.numeric(intVal_ab[2])*1.18, as.numeric(intVal_ab[2])*1.18, length =(varyingYpoint_mid2 - varyingYpoint_mid - 1)),
      seq(as.numeric(intVal_ab[2])*1.18, as.numeric(intVal_ab[2])*1.18*1.05, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(intVal_ab[2])*1.18*1.05, POC_AU$npts - varyingYpoint_end + 1))
  
  dfList$tau_ab[5, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_ab[3]),as.numeric(intVal_ab[3]) , length = (varyingYpoint_fir - 2)),
      seq(as.numeric(intVal_ab[3]), as.numeric(intVal_ab[3])*1.18, length =(varyingYpoint_mid - varyingYpoint_fir + 1)),
      rep(as.numeric(intVal_ab[3])*1.18, as.numeric(intVal_ab[3])*1.18, length =(varyingYpoint_mid2 - varyingYpoint_mid - 1)),
      seq(as.numeric(intVal_ab[3])*1.18, as.numeric(intVal_ab[3])*1.18*1.05, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(intVal_ab[3])*1.18*1.05, POC_AU$npts - varyingYpoint_end + 1))
  
  dfList$tau_RNA[3, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_RNA[1]),as.numeric(intVal_RNA[1]) , length = (varyingYpoint_fir - 2)),
      seq(as.numeric(intVal_RNA[1]), as.numeric(intVal_RNA[1])*1.08, length =(varyingYpoint_mid - varyingYpoint_fir + 1)),
      rep(as.numeric(intVal_RNA[1])*1.08, as.numeric(intVal_RNA[1])*1.08, length =(varyingYpoint_mid2 - varyingYpoint_mid - 1)),
      seq(as.numeric(intVal_RNA[1])*1.08, as.numeric(intVal_RNA[1])*1.08*1.0056, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(intVal_RNA[1])*1.08*1.0056, POC_AU$npts - varyingYpoint_end + 1))
  
  dfList$tau_RNA[4, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_RNA[2]),as.numeric(intVal_RNA[2]) , length = (varyingYpoint_fir - 2)),
      seq(as.numeric(intVal_RNA[2]), as.numeric(intVal_RNA[2])*1.08, length =(varyingYpoint_mid - varyingYpoint_fir + 1)),
      rep(as.numeric(intVal_RNA[2])*1.08, as.numeric(intVal_RNA[2])*1.08, length =(varyingYpoint_mid2 - varyingYpoint_mid - 1)),
      seq(as.numeric(intVal_RNA[2])*1.08, as.numeric(intVal_RNA[2])*1.08*1.0056, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(intVal_RNA[2])*1.08*1.0056, POC_AU$npts - varyingYpoint_end + 1))
  
  dfList$tau_RNA[5, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal_RNA[3]),as.numeric(intVal_RNA[3]) , length = (varyingYpoint_fir - 2)),
      seq(as.numeric(intVal_RNA[3]), as.numeric(intVal_RNA[3])*1.08, length =(varyingYpoint_mid - varyingYpoint_fir + 1)),
      rep(as.numeric(intVal_RNA[3])*1.08, as.numeric(intVal_RNA[3])*1.08, length =(varyingYpoint_mid2 - varyingYpoint_mid - 1)),
      seq(as.numeric(intVal_RNA[3])*1.08, as.numeric(intVal_RNA[3])*1.08*1.0056, length =(varyingYpoint_end - varyingYpoint_mid2 + 1)),
      rep(as.numeric(intVal_RNA[3])*1.08*1.0056, POC_AU$npts - varyingYpoint_end + 1))
  
}


varying_Yint <- 2015
varying_Yfir <- 2018
varying_Ymid <- 2019
varying_Yend <- 2021
calibration_Y <- 2015
varyingYpoint_int <- (varying_Yint - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_fir <- (varying_Yfir - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_mid <- (varying_Ymid - calibration_Y)/POC_AU$timestep + 1
varyingYpoint_end <- (varying_Yend - calibration_Y)/POC_AU$timestep + 1

for ( i in 2:dim(dfList$eta)[[2]]){  
  dfList$eta[3, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[3]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, 0.9995, length = (varyingYpoint_mid - varyingYpoint_fir)),
      rep(0.9995, POC_AU$npts - varyingYpoint_mid))
  
  dfList$eta[4, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[3]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, 0.9995, length = (varyingYpoint_mid - varyingYpoint_fir)),
      rep(0.9995, POC_AU$npts - varyingYpoint_mid))
  
  
  
  dfList$eta[5, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[3]), 0.6, length = (varyingYpoint_fir)),
      seq(0.6, 0.9995, length = (varyingYpoint_mid - varyingYpoint_fir)),
      rep(0.9995, POC_AU$npts - varyingYpoint_mid))
  
  
  dfList$rho[1, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(0, 0.5, length = (varyingYpoint_fir)),
      seq(0.5, 0.5, length = (varyingYpoint_mid - varyingYpoint_fir)),
      rep(0.5, POC_AU$npts - varyingYpoint_mid))
  
  dfList$rho[2, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(0, 0.5, length = (varyingYpoint_fir)),
      seq(0.5, 0.5, length = (varyingYpoint_mid - varyingYpoint_fir)),
      rep(0.5, POC_AU$npts - varyingYpoint_mid))
  
  dfList$rho[3, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(0, 0.6, length = (varyingYpoint_fir)),
      seq(0.6, 0.999, length = (varyingYpoint_mid - varyingYpoint_fir)),
      rep(0.999, POC_AU$npts - varyingYpoint_mid))
  
  
  dfList$rho[4, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(0, 0.6, length = (varyingYpoint_fir)),
      seq(0.6, 0.999, length = (varyingYpoint_mid - varyingYpoint_fir)),
      rep(0.999, POC_AU$npts - varyingYpoint_mid))
  
  
  
  
  dfList$rho[5, i, c(varyingYpoint_int:POC_AU$npts)] <- 
    c(seq(0, 0.6, length = (varyingYpoint_fir)),
      seq(0.6,0.999, length = (varyingYpoint_mid - varyingYpoint_fir)),
      rep(0.999, POC_AU$npts - varyingYpoint_mid))
  
  
}

intVal <- dfList$cured[, 3, 1] 

varying_Yfir <- 2016

calibration_Y <- 2015

varyingYpoint_fir <- 2*(varying_Yfir - calibration_Y)/POC_AU$timestep


for ( i in 1:dim(dfList$cured)[[1]]){  
  dfList$cured[i, "f0", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
  
  dfList$cured[i, "f1", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
  
  dfList$cured[i,"f2", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
  
  dfList$cured[i, "f3", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
  
  dfList$cured[i,"f4", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
  
  dfList$cured[i,"dc", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
  dfList$cured[i,"hcc", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
  dfList$cured[i,"lt", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
  dfList$cured[i,"plt", c(1:POC_AU$npts)] <- 
    c(seq(as.numeric(intVal[i]), 0.9868763, length = (varyingYpoint_fir)),
      rep(0.9868763, POC_AU$npts - varyingYpoint_fir))
}
# treatment duration 
dfList$lota <- dfList$cured


tic <- proc.time()
endY <- 100
calibrateInit <- HCVMSM(POC_AU, best_estimates, best_est_pop,
                        disease_progress,pop_array,
                        dfList,dfList_sc ,fib, 
                        modelrun="UN", proj = "POC_AU", end_Y = endY, 
                        fc_sc = fc,
                        fp = NULL)


toc <- proc.time() - tic 
toc



save(calibrateInit, bvalue, TreatInit,best_estimates, dfList,
     file = file.path(OutputFolder ,
                      paste0(project_name,"cali_timev" ,".rda")))

# quick check out the key indicators 
endY_plot <- 2030- POC_AU$cabY

subpop_N <- lapply(POC_AU$popNames, function(x){ 
  
  a <- popResults_MidYear(POC_AU, calibrateInit,
                          Population = x,
                          Disease_prog = NULL, 
                          Cascade = NULL, param = NULL, 
                          endYear = endY)%>%ungroup()
})

names(subpop_N) <- POC_AU$popNames
pop_N <- dplyr::bind_rows(subpop_N, .id = 'population')
pop_N <- pop_N%>%arrange(year, population)
tempNOTInfectedRNA_subpop <- popResults_MidYear(POC_AU, calibrateInit,Population = POC_AU$popNames,
                                                Disease_prog = NULL, 
                                                Cascade = c("s", "cured"), 
                                                param = NULL ,endYear = endY)%>%
  as_tibble()%>%group_by(year, population)%>%
  summarise(best = sum(best))%>%arrange(year, population)



tempPrevRNA_subpop <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), each = POC_AU$npops),
                            population = POC_AU$popNames,
                            
                            as.data.frame(100*(pop_N[, -c(1,2)] - 
                                                 tempNOTInfectedRNA_subpop[ ,-c(1,2)])/ 
                                            pop_N[ ,-c(1,2)]))%>%
  tibble::as_tibble()  
pop_labname <- c("PWID in community",  "Former PWID in community", 
                 "PWID in prisons",  "Former PWID in prisons", 
                 "nonPWID in prisons")
tempPrevRNA_subpop <- tempPrevRNA_subpop%>%
  mutate(population = factor(population, 
                             levels = POC_AU$popNames, 
                             labels = pop_labname ))
HCVPrevRNA <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVPrevRNA_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = HCV.RNA.prevalence*100,
                           up = upper*100,
                           low = lower*100,
                           population = factor(population, 
                                               levels = POC_AU$popNames, 
                                               labels = pop_labname ))
popPrevRNAPlot <- indicatorPlot(POC_AU, tempPrevRNA_subpop, 
                                ylabel = "HCV prevalence (%)",
                                xlimits = c(POC_AU$startYear, 
                                            (POC_AU$startYear+endY_plot), 5),
                                calibration_Y = POC_AU$cabY,
                                rangeun = NULL, 
                                groupPlot = NULL, 
                                facetPlot = population,
                                observationData = HCVPrevRNA , 
                                simulateYear = POC_AU$simY) +
  ggtitle("HCV RNA prevalence by population") 

popPrevRNAPlot <- popPrevRNAPlot + 
  facet_custom (~population,
                scales = "free", ncol = 3,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 60))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 60))),
                    
                    scale_new(3,
                              scale_y_continuous(limits = 
                                                   c(0, 80))),
                    scale_new(4,
                              scale_y_continuous(limits = 
                                                   c(0, 80))),
                    scale_new(5,
                              scale_y_continuous(limits = 
                                                   c(0, 20)))
                  )) + theme_bw()


calibrateFlow <- calibrateInit[!names(calibrateInit)%in%
                                 c("allPops", "newpop_tran", "newpop_tranState", "HCVdeathState",
                                   "newDeathState", "death_hcv")]
flow_sub <- list()


flow_sub <- lapply(names(calibrateFlow), function(x){ 
  a <- indicatorResults(POC_AU, calibrateFlow, x, 
                        pop=POC_AU$popNames,
                        paramR = NULL, range = NULL,
                        endY = endY)
})

names(flow_sub) <- names(calibrateFlow)


flow_setting <- lapply(flow_sub, function(x){ 
  
  a <- x%>%
    mutate(setting = ifelse(population %in% c("C_PWID", "C_fPWID"), 
                            "commu", "prisons"))
  
  a <- a%>%group_by(year, setting)%>%summarise_at("best", sum)%>%
    mutate(population = setting)%>%select(-setting)
}) 
N_treatment <- cbind(year = rep(seq(POC_AU$startYear , endY-1 ,1), 
                                each = 2),
                     population = flow_setting$newTreatment[ ,3],
                     as.data.frame(flow_setting$newTreatment[, -c(1,3)] + 
                                     flow_setting$newRetreat[, -c(1,3)] + 
                                     flow_setting$newTreatment_sc[, -c(1,3)]))%>%
  tibble::as_tibble()%>%
  mutate(population = factor(population, 
                             levels = c("commu", "prisons"), 
                             labels = c("Community", "Prisons" )))



HCVtreatinitN_setting_fit <-read.csv(file.path(paste0(DataFolder%>%dirname(), "/HCVtreatinitN_setting_POC_AU.csv")), header = TRUE)%>%
  as.data.frame()%>%mutate(time = year- POC_AU$cabY + 1, 
                           realPop = realpop,
                           up = upper,
                           low = lower,
                           population = factor(population, 
                                               levels = c("commu", "prisons"), 
                                               labels = c("Community", "Prisons")))


N_treatment_setting_p <- indicatorPlot(POC_AU, N_treatment, 
                                       ylabel = "N",
                                       xlimits = c(POC_AU$startYear, 
                                                   (POC_AU$startYear+endY_plot), 5),
                                       calibration_Y = POC_AU$cabY,
                                       rangeun = NULL, 
                                       groupPlot = NULL, 
                                       facetPlot = population,
                                       observationData = HCVtreatinitN_setting_fit, 
                                       simulateYear = POC_AU$simY) + 
  ggtitle("Number of treatment init") + theme_bw()

N_treatment_setting_p <- N_treatment_setting_p + 
  facet_custom (~population,
                scales = "free", ncol = 2,
                scale_overrides = 
                  list(
                    scale_new(1,
                              scale_y_continuous(limits = 
                                                   c(0, 40000))),
                    scale_new(2,
                              scale_y_continuous(limits = 
                                                   c(0, 5000)))))
popPrevRNAPlot

N_treatment_setting_p
tempPrevRNA_subpop%>%filter(year%in% c(8,9,16))
N_treatment%>%filter(year %in% c(7,8,9,10))
